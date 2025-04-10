C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C


      subroutine treestore(numpom)
#include "aaa.h"
      parameter (mx3ptl=150000)
      integer num (mxptl+mx3ptl)
      parameter (mxptlxx=1.75*mxptl)
      integer id (mxptlxx), ist (mxptlxx), ity (mxptlxx)
      integer ior (mxptlxx), jor (mxptlxx)
      real px (mxptlxx), py (mxptlxx), pz (mxptlxx), eam (mxptlxx)
      real x (mxptlxx), y (mxptlxx), z (mxptlxx), t (mxptlxx)
      real zus (mxptlxx)
      integer iversnx,laprojx,maprojx,latargx,matargx,nfullx,nfreezex
      real engyx, sigt, dy
      integer nhard, npartproj, nparttarg 
      integer nspecprojp, nspecprojn, nspectargp, nspectargn
      common/cNpart/sumNpart
      common/ciniflo/ecccoe(5),eccphi(5)
      common/cigrpac/igrpac/cigrval/igrvalF
      common/cnfrx/nfrx
      common/cmxpom/mxpom
      character*200 fn

      if(irootcproot.gt.0.or.hepmc.eq.1)write(ifmt,'(a)')
     *'enter treestore'

      itest=0
 
      call maxsize_get(2, mx2ptl )
      if(mx3ptl.ne.mx2ptl)stop'ERROR 07052022'

      nfillt(numpom)=nfillt(numpom)+1

      if(nfillt(numpom).eq.1)then
        nopent(numpom)=nopent(numpom)+1
        if(mxpom.eq.0)then !------normal case-----
          if(numpom.ne.1)stop'ERROR 18042014'
          ii=index(fndt(1:nfndt),".data")-1
          if(irootcproot.gt.0)then
            write(ifmt,'(a)')'opentree '//fndt(1:ii)//'.root'
            call opentree(1,iextree,ihepmc3,
     .                    fndt(1:ii)//'.root'//CHAR(0))
          endif
          if(hepmc.eq.1)then
            write(ifmt,'(a)')'openhepmc '//fnhm(1:nfnhm)//'.hepmc'
            call openhepmc(fnhm(1:nfnhm)//'.hepmc'//CHAR(0))
          else if(hepmc.eq.2)then
            write(ifmt,'(a)')'openhepmc /dev/stdout'
            call openhepmc('/dev/stdout'//CHAR(0))
          endif 
          izm=ifillTree
          nfullt(1)=nfull
        else !-----orderTree case using Poms or similar-----
          call tfname(numpom,fn,nfn)          
          if(irootcproot.gt.0)then
            write(ifmt,'(a,2i7)')'opentree '//fn(1:nfn)
     .     ,numpom,nfullt(numpom)
            call opentree(numpom,iextree,ihepmc3,fn(1:nfn)//CHAR(0))
          endif
          if(hepmc.eq.1)then
            call openhepmc(fnhm(1:nfnhm)//'.hepmc'//CHAR(0))
          else if(hepmc.eq.2)then
            call openhepmc('/dev/stdout'//CHAR(0))
          endif
          izm=izmode
        endif
        iversnx=iversn
        laprojx=laproj
        maprojx=maproj
        latargx=latarg
        matargx=matarg
        engyx  =engy
        nfullx =nfullt(numpom)
        nfreezex=nfreeze
        if(irootcproot.gt.0)
     .  call FillHead(numpom,iversnx,laprojx,maprojx,latargx,matargx
     .  ,engyx,nfullx,nfreezex)
        sumNpart=0
      endif

      call maxsize_get(2, mx2ptl ) 
      do i=1,mxptl+mx2ptl
        num(i)=-1
      enddo

      sumNpart=sumNpart+nptl
      n=0
      
      do iloo=1,nptl
      
        i=iloo   
        if(iloo.gt.mxptl)then
          call restorecccptl(iloo,mxptl+2)
          i=mxptl+2
        endif

        jo=jorptl(i)
        !if(istptl(i).eq.20.or.istptl(i).eq.21)jo=nint(radptl(i)) !KW2309 makes problems
        !print*,'          ',iorptl(i),jo,i,'   '
        !.   ,idptl(i)
        if(istptl(i).le.1.or.istptl(i).eq.29
     .   .or.istptl(i).eq.21.or.istptl(i).eq.25
     .   .or.istptl(i).eq.26
     .   .or.istptl(i).eq.3.or.istptl(i).eq.8
     .   .or.istptl(i).eq.9.or.istptl(i).eq.6
     .   .or.istptl(i).eq.-1.or.istptl(i).eq.-2
     .   .or.istptl(i).eq.31
     .   .or.istptl(i).eq.7)then
          n=n+1
          if(n.gt.mxptlxx)
     .    stop'####### ERROR : mxptlxx too small #######'
          num(iloo)=n
          id(n)= idptl(i)
          ist(n)=istptl(i)
          ity(n)=ityptl(i)
          io=iorptl(i)
          if(istptl(i).le.1)then
            if(io.gt.0)then
              zus(n)=num(io)      !mother identified (-1 if not stored)
            elseif(io.eq.-999)then
              zus(n)=-999          !decay product from decay inside cascade
            else
              zus(n)=-2    !mnu=-2 -> no mother
            endif
          elseif(istptl(i).eq.25)then
            zus(n)=rinptl(i)             ! DelE
            ity(n)=nint(qsqptl(i)*1000)  ! L*1000
          else
            zus(n)=0  !presently unused
          endif
          if(io.gt.0)then
            ior(n)=max(0,num(io))
          else
            ior(n)=0
          endif
          if(jo.gt.0)then
            jor(n)=max(0,num(jo))
          else
            jor(n)=0
          endif
          !if(ist(n).eq.21.and.jor(n).gt.0)ior(jor(n))=ior(n) !KW2309 makes problems
          px(n)= pptl(1,i)
          py(n)= pptl(2,i)
          pz(n)= pptl(3,i)
          eam(n)= pptl(5,i) !pptl(4,i) !changed in 3.127
          x(n)= xorptl(1,i)
          y(n)= xorptl(2,i)
          z(n)= xorptl(3,i)
          t(n)= xorptl(4,i)
        endif
        
      enddo
      
      if(maproj+matarg.gt.10)then
        !use zus(n) to store additional event information
        !with no need to change tree structure
        do m=2,5
          zus(  m-1)=eccphi(m)
          zus(4+m-1)=ecccoe(m)
        enddo
        !if(nfrx.eq.0)
        !.  write(ifmt,'(a,4f6.2,3x,4f7.4)')' ++treestore++ eccphi/coe='
        !.  ,(eccphi(m),m=2,5),(ecccoe(m),m=2,5)
      endif
      nptevt=n
      iorpre=0
      !do n=1,nptevt
      !  if(ist(n).eq.21)then
      !    if(ior(n).ne.iorpre)then
      !      sumx=0
      !      sumy=0
      !    endif
      !    sumx=sumx+px(n)
      !    sumy=sumy+py(n)
      !    iorpre=ior(n)
      !  else
      !    sumx=0
      !    sumy=0
      !  endif
      !  print*,'treestore ',ior(n),jor(n),n,'***'
      !.    ,id(n),ist(n),ity(n),'   '
      !.    ,sqrt(px(n)**2+py(n)**2)
      !.    ,sqrt(sumx**2+sumy**2)
      !enddo
      !---------------------------------------------------------
      !write(ifch,*)'==============TreeWrite==========='
      !do n=1,nptevt
      !if(ist(n).eq.0)write(ifch,*)n,id(n),px(n),py(n),pz(n)
      !enddo
      !---------------------------------------------------------
      call getevt(nev,phi,phir,bim,egy,npt,ngl,kol)
      call getecc(psi2,psi3,psi4,psi5,ecci2,ecci3,ecci4,ecci5)
      cnt=centrVar(izm)
      sigt=sigtot
c*JJ      nhard=ikoevt
      nhard=0 !*JJ ikoevt is number of pomerons, so not exactly nhard
      npartproj=ng11evt !*JJ number of Glauber proj participants
      nparttarg=ng12evt !*JJ number of Glauber targ participants
c     nspecp=ngspecp !nppevt+ntpevt
c     nspecn=ngspecn !npnevt+ntnevt
c     nspecp=nppevt+ntpevt
c     nspecn=npnevt+ntnevt
      nspecprojp=nppevt
      nspecprojn=npnevt
      nspectargp=ntpevt
      nspectargn=ntnevt
      if(irootcproot.gt.0)then
        call fillTree(numpom,nptevt, cnt, sigt, iextree, ihepmc3
     .  , nev, npt, ngl, kol, nhard, npartproj, nparttarg
     .  , nspecprojp, nspecprojn, nspectargp, nspectargn 
     .  , phi, phir, psi2, psi3, psi4, psi5, ecci2, ecci3, ecci4, ecci5
     .  , id, ist, ity, ior, jor, zus, px, py, pz, eam 
     .  , x, y, z, t, iret)
      endif 
      if(hepmc.eq.1 .or. hepmc.eq.2)then
        if(ish.ge.2)call alist('Fill HepMC direct&',0,0)
        call getReaction(iprojZ,iprojA,itargZ,itargA,fegyevt)
        dy=-hepmc_rapcms
        call fillhepmc(iextree,ihepmc3,nrevt,fegyevt,dy, iprojZ, iprojA
     .  , itargZ, itargA, nhard, ngl, npartproj, nparttarg
     .  , nspecprojp, nspecprojn, nspectargp, nspectargn, nptevt, cnt
     .  , sigt, id, ist, ity, ior, px, py, pz, eam, x, y, z, t
     .  , hepmc_record_mode, hepmc_tau_decay, hepmc_record_id_nb
     .  , hepmc_record_id_list)

      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !  List of fillTree variables
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! nptevt ....... number of particles in event
      ! cnt .......... centrality variable
      ! sigt ......... total cross-section of the process
      ! iextree ...... tree extension (0,1)
      ! ihepmc3 ...... HepMC format version (0=HepMC2,1=HepMC3)
      ! id ........... particle id (epos code, 
      !                  see "List of particle codes" in "KW/ids.f")
      ! ist .......... status and particle type
      !                 0 and 1 ... hadrons (0 = last generation)
      !                 3 ......... intermediate, used when hacas is called
      !             As a "special service", in case of a full simulation
      !             (namely EPOS with hydro and hadronic cascade (hacas)) 
      !             in addition to the list of particles with ist 0 and 1
      !             we provide another list, which corresponds to EPOS
      !             with hydro but without hacas, with the following ist:         
      !                 8 and 6 ... without hacas (8 = last generation), 
      !                               used when hacas is called 
      !             Since some procedures cannot guarantee strict energy 
      !             conservation, we provide an additional list of particles 
      !             where a correction procedure has been employed to 
      !             enforce energy conservation, with the following ist:
      !                -2 and -1 .. Energy strictly conserved  
      !                             (-2 = last generation)
      !                 21 ........ partons
      !                 25 ........ intermediate out-Born partons
      !                 29 ........ string
      ! ity .......... type of particle origin
      !                 20-29 Soft Pomeron
      !                 30-39 Semihard Pomeron
      !                 40-49 Projectile remnants
      !                 50-59 Target remnants
      !                 60,61 Plasma
      ! ior, jor ..... indices of parents
      !                 ior>0, jor=0 : ior is the mother index
      !                 ior>0, jor>0 : parent indices from  ior to jor
      !                        (a string has many parents (partons))
      !                 ior=0, jor=0 : no parent known
      !                 With hacas, the parent information is lost, but there 
      !                 are special lists of resonances for particular 
      !                 decay channels. Contact authors.   
      ! px, py, pz ... particle  momentum 
      ! eam .......... mass (before 3.127 energy)
      ! x, y, z, t ...... particle four position
      ! zus ............. for private use 
      !                   (the meaning may change from version to version)
      !                   in case of partons : presently unused
      !                   in case of hadrons :  decay information :
      !                    -999 : hadron is decay product from decay
      !                           in cascade part (mother unknown)
      !                      -1 : hadron is decay product, mother not stored
      !                      >0 : hadron is decay product, mother index = zus
      !                      -2 : no mother
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(iret.ne.0)
     .stop'##### ERROR in treestore: array in fillTree too small #####'
      if(irootcproot.gt.0.or.hepmc.eq.1)write(ifmt,'(a)')
     *'exit treestore'
      end
      subroutine tfname(numpom,fn,nfn)          
#include "aaa.h"
      character cfmt*4,ctxt*8
      common/cigrpac/igrpac/cigrval/igrval
      common/cmxpom/mxpom
      character*200 fn
      cfmt='(i_)'
      ii=index(fndt(1:nfndt),".data")-1
      if(mxpom.eq.0)then
        nfn=ii+5
        fn=fndt(1:ii)//'.root'
      else  
        do j=ii-1,1,-1
         if(fndt(j:j).eq.'-')goto1
        enddo
   1    jj=j
        read(fndt(jj+1:ii),*)kkx
        kk=(kkx-1)*igrpac+igrval
        nn=(kk-1)*mxpom+numpom
        ndig=1+int(log10(float(nn)))
        write(cfmt(3:3),'(i1)')ndig
        write(ctxt(1:ndig),cfmt)nn
        nfn=jj+ndig+5
        fn=fndt(1:jj)//ctxt(1:ndig)//'.root'
      endif
      end

      !-----------------------------------------------------------------
      subroutine pthardparent(i,i1,i2,ptmax)
#include "aaa.h"
      i1=0
      i2=0
      ptmax=0
      if(istptl(i).gt.1)return
      io=i
      iso=0
      if(io.gt.0)iso=istptl(io)
      do while (io.gt.0.and.iso.ne.29)
        io=iorptl(io)
        iso=0
        if(io.gt.0)iso=istptl(io)
      enddo
      if(io.le.0)then
        !write(ifch,*)' nonpositive io  ', i,io
        return
      elseif(istptl(io).eq.29)then
         ptmax=0
         do k=iorptl(io),jorptl(io)
           ih=nint(radptl(k))
           pt=0
           if(ih.gt.0)then
             if(i1.eq.0)then
               i1=ih
             elseif(i2.eq.0.and.ih.ne.i1)then
               i2=ih
             endif
             pt = sqrt ( pptl(1,ih)**2 + pptl(2,ih)**2 )
             ptmax=max(ptmax,pt)
           endif
           !write(ifch,*)' string  ',k,ih,pt,ptmax
         enddo
      else
         stop'\n\n STOP in pthardparent\n\n'
      endif
      end

      subroutine treeclose
      common /cnnnhis/nnnhis
#include "aaa.h"
      character*200 fn
      do npom=1,mxxpom
       if(nopent(npom).gt.0)then
        if(irootcproot.gt.0)then
          write(ifmt,'(a,i7)')'closetree',npom
          call clop(3)
          call tfname(npom,fn,nfn)          
          call closetree(npom,iextree,ihepmc3,fn(1:nfn)//CHAR(0))
          write(ifmt,'(a,i7)')'closetree done'
          call clop(3)
        endif
        if(hepmc.eq.1 .or. hepmc.eq.2)then
          write(ifmt,'(a)')'closehepmc'
          call clop(3)
          call closehepmc()
          write(ifmt,'(a,i7)')'closehepmc done'
          call clop(3)
        endif
       endif
      enddo
#if !__BS__ 
      !if(ifemto.gt.0)call femto !femto removed in 3444m4
      if(ivmd.gt.0)call vmd     !maybe better move to epos-bas (like for dih)
#endif
      if(igrTree.gt.0.and.nnnhis.gt.0
     ..or.irootcproot.eq.01)call getCevtCptl
      end

c###################################################################################
c###################################################################################

      !----------------------------
      subroutine getNames(mode)
      !----------------------------
#include "aaa.h"
      common/cNpart/sumNpart
      !                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*400       fnamemtr,    cbasout
      common/croot/imtr,fnamemtr /croot6/cbasout,iba,nout
      !                 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      common/croot3/iMaxNpart
      common/croot5/jj1,jj2
      character cext3*10 
      common/ccext3/iext3,cext3 
      imtr=ifmt
      ii=index(fnhi(1:nfnhi),".histo")-1
      do j=ii-1,1,-1
       if(fnhi(j:j).eq.'-')goto1
      enddo
  1   jj=j
      if(jj.eq.1)then
        write(ifmt,'(//a/a/a,a/a,i10,i10)')
     .   'STOP in getFileNames'
     .  ,'problem with histo filename'
     .  ,'filename = ',fnhi(1:nfnhi)
     .  ,'jj,ii = ',jj,ii
        stop
      endif
      jj=jj+1
      !write(ifmt,*)'+++++',fnhi(1:nfnhi),jj,ii
      nout=0
      cbasout=fnhi(1:ii)//' '
      if(jj.gt.2)then
        read(fnhi(jj:ii),*)nout
        cbasout=fnhi(1:jj-2)//' '
      endif
      iba=index(cbasout,' ')-1

      if(mode.eq.3)return

      jj2=jj-2
      if(cext3(1:iext3).ne.'-')jj2=jj2-iext3
      do j=jj2,1,-1
       if(fnhi(j:j).eq.'/')goto2
      enddo
  2   jj1=j+1
      fnamemtr=fnmt(1:nfnmt)//' '
      if(ifmt.eq.6)fnamemtr=' '

      if(mode.eq.1)return

      if(igrTree.gt.0)then
      iMaxNpart=muTree
      else
      iMaxNpart=sumNpart/nfull/nfreeze
      endif

      end

      !----------------------------
      subroutine checkRootFiles(ii,i1,i3,i4)
      !----------------------------
#include "aaa.h"
      logical rootfile(10000)
      common/croot4/rootfile
      !              ~~~~~~~~~~~~~~~
      character*400 fnamein
      common/croot2/fnamein,ix1
      !              ~~~~~~~~~~~~~~~
      common/croot5/jj1,jj2
      common/ciext4/iext4
      character*4 cfmt, cij
      character*1 c4
      c4='0'
      ic4=iext4
      rootfile(1)=.false.
      if(ifillTree.gt.0)rootfile(1)=.true.
      imax=max(1,igrTree)
      i1=1
      i2=imax
      i3=1
      i4=imax
      if(ifillTree.le.0.and.igrTree.le.0)then
        write(ifmt,'(a//a//)')'STOP','ifillTree <= 0  & igrTree <= 0'
        stop
      endif
      if(igrTree.gt.0)then
        if(iopcnt.le.0)stop'\n\n batch required \n\n'
        i1=(iopcnt-1)*imax+1
        i2=iopcnt*imax
        cfmt='(i )'
        i3=0
        i4=-1
        if(ii.eq.2)then
        write(ifmt,'(a/7x,a/7x,a)')'checkRootFiles: no valid file found'
     .   ,'in the following list, F means missing file, '
     .   ,'  T means other problem, see error file'
        else   
        write(ifmt,'(2(a,i8))')'checkRootFiles: loop',i1,'  to',i2
        endif
        do i=i1,i2
         ij=(i-1)/10000+1
         nij=1+log10(float(ij))
         if(ij.le.9)then
          write(cij(1:nij),'(i1)')ij
         elseif(ij.le.99)then
          write(cij(1:nij),'(i2)')ij
         elseif(ij.le.999)then
          write(cij(1:nij),'(i3)')ij
         elseif(ij.le.9999)then
          write(cij(1:nij),'(i4)')ij
         else
          stop'\n\n ERROR 27042011\n\n'
         endif
         if(isyst.eq.0)then       !----system i (interactive local)---
           fnamein=
c    .     fnnx(1:nfnnx)//'../../../public/root/'
     .     '/sps/nantheo/werner/ss/' !*JJ
     .     //fnhi(jj1+2:jj2)//c4(1:ic4)
     .     //"/"//cij(1:nij)//"/"//fnhi(jj1:jj2)//c4(1:ic4)//' '
         elseif(isyst.eq.1)then       !-------system q (qepos)---------
           fnamein='/tmp/'//fnhi(jj1:jj2)//c4(1:ic4)//' '
         elseif(isyst.eq.2)then       !-------system c (cepos)---------
           fnamein='./'//fnhi(jj1:jj2)//c4(1:ic4)//' '
         elseif(isyst.eq.-1)then     !--system -q (interactive q)
           ju=2
           fnamein=
     .     '/scratch/theoric2/werner/public/root/'
     .     //fnhi(jj1+ju:jj2)//c4(1:ic4)
     .     //"/"//cij(1:nij)//"/"//"/z-"//fnhi(jj1+ju:jj2)
     .     //c4(1:ic4)//' '
         elseif(isyst.eq.-2)then      !--system -c (interactive c)
           ju=2
           fnamein=
     .     '/sps/nantheo/werner/ss/'//fnhi(jj1+ju:jj2)//c4(1:ic4)
     .     //"/"//cij(1:nij)//"/"//"/z-"//fnhi(jj1+ju:jj2)
     .     //c4(1:ic4)//' '
         endif
         ix1=index(fnamein,' ')
         fnamein(ix1:ix1)='-'
         n=1+log10(float(i))
         write(cfmt(3:3),'(i1)')n
         write(fnamein(ix1+1:ix1+n),cfmt)i
         fnamein(ix1+n+1:ix1+n+6)='.root '
         inquire(file=fnamein(1:index(fnamein,' ')-1)
     .         ,exist=rootfile(i-i1+1))
         if(ii.eq.2)then
           write(ifmt,*)i,'  ',fnamein(1:index(fnamein,' ')-1),'  '
     .     ,rootfile(i-i1+1)
         endif
         if(rootfile(i-i1+1))then
          call checktreefile(
     .    fnamein(1:index(fnamein,' ')-1)//CHAR(0),ierr)
          if(ierr.ne.0)rootfile(i-i1+1)=.false.
          if(rootfile(i-i1+1).and.i3.eq.0)i3=i
          if(rootfile(i-i1+1))i4=i
         endif
        enddo
      endif
      end

      !----------------------------
      subroutine updateName(i)
      !----------------------------
#include "aaa.h"
      !              ~~~~~~~~~~~~~~~
      character*400  fnamein
      character*400  fnameinh
      common/croot2/fnamein,ix1
      !              ~~~~~~~~~~~~~~~
      character*4 cfmt
      if(igrTree.gt.0)then
        n=1+log10(float(i))
        cfmt='(i )'
        write(cfmt(3:3),'(i1)')n
        write(fnamein(ix1+1:ix1+n),cfmt)i
        fnamein(ix1+n+1:ix1+n+12)='.root       '
        fnameinh(ix1+n+1:ix1+n+12)='.hepmc      '
      endif
      end

c###################################################################################
      subroutine getCevtCptl
#include "aaa.h"
      character*400 fnamemtr,fnamein, cbasout
      logical rootfile(10000)
      common/croot/imtr,fnamemtr /croot2/fnamein,ix1 
      common/croot6/cbasout,iba,nout
      common/croot4/rootfile
      common /cnnnhis/nnnhis
      common/cirtfile/irtfile
      parameter (mxptlxx=1.75*mxptl)
      integer id (mxptlxx), ist (mxptlxx), ity (mxptlxx)
      integer ior (mxptlxx), jor (mxptlxx)
      real px (mxptlxx), py (mxptlxx), pz (mxptlxx), eam (mxptlxx)
      real p5 (mxptlxx), en (mxptlxx)
      real x (mxptlxx), y (mxptlxx), z (mxptlxx), t (mxptlxx)
      real sigt, dy
      integer nhard, npartp, npartt 
      integer nspecprojp, nspecprojn, nspectargp, nspectargn
      data nrevtxx/0/
      save nrevtxx

      write(ifmt,'(67a)')('#',i=1,26),' getCevtCptl ',('#',i=1,28)
      call clop(3)

      call getNames(2)
      call checkRootFiles(1,i1,i3,i4)
      if(i4.eq.-1)call checkRootFiles(2,i1,i3,i4)
      write(ifmt,'(2(a,i8))')'getCevtCptl: loop',i3,'  to',i4
      if(hepmc.eq.1)then
        call openhepmc(fnhm(1:(nfnhm))//'.hepmc'//CHAR(0))
      else if(hepmc.eq.2)then
        call openhepmc('/dev/stdout'//CHAR(0))
      endif
      do i=i3,i4
        call clop(3)
        call updateName(i)
        if(rootfile(i-i1+1))then
         call etreeopen(fnamein(1:index(fnamein,' ')-1)//CHAR(0))
         call etreenevent(nevents)
         call etreehead(iversn,laproj,maproj,latarg,matarg,engy
     .   ,nfreeze,nfull)
         egyevt=engy
#if !__BS__ 
         call transfer_head(iversn,laproj,maproj,latarg,matarg,engy)
#endif
         do ii=1,nevents
           nrevtxx=nrevtxx+1
           nrevt=nrevtxx
           call etreeevt(ii,np,cnt,sigt,iextree,ihepmc3
     .     ,nev,npt,ngl,kol,nhard,npartp,npartt,nspecprojp,nspecprojn
     .     ,nspectargp,nspectargn,phi,phir,psi2,psi3,psi4,psi5
     .     ,ecci2,ecci3,ecci4,ecci5,id,ist,ity,ior,jor,px,py,pz,eam
     .     ,x,y,z,t)
           nptl=np
           ng1evt=npartp+npartt
           nglevt=ngl
           irtfile=i
           call setcentrVar(cnt)
           !write(ifmt,*)'getCevtCptl np = ',np
           !print*,'+++++',nptl,bimevt,izmode,i
           iold=0
           if(nptlpt.eq.maproj+matarg)then
             if(abs(px(1)).lt.1e-5.and.abs(py(1)).lt.1e-5
     .       .and.sqrt(pz(1)**2+0.94**2).gt.0.50*engy/2)then
               if(eam(1).gt.0.50*engy/2)iold=1 !looks like an energy and not a mass
             else
               if(iphsd.eq.0)then
                 print*,id(1),ist(1),ity(1),px(1),py(1),pz(1),eam(1) 
                 stop'####### ERROR 26032016 #######'
               endif
             endif
           endif
           do n=1,np
             ior(n) = ior(n) + 1
             jor(n) = jor(n) + 1
             call checkcccptl(n)
             nxx=n
             if(n.gt.mxptl)nxx=mxptl+1
             idptl(nxx)    = id(n)
             istptl(nxx)   = ist(n)
             ityptl(nxx)   = ity(n)
             iorptl(nxx)   = ior(n)
             jorptl(nxx)   = jor(n)
             pptl(1,nxx)   = px(n)
             pptl(2,nxx)   = py(n)
             pptl(3,nxx)   = pz(n)
             if(iold.eq.1)then
               e = eam(n)
               if(ity(n).eq.61)then
                 stop'ERROR 25032021'  !call getQBM(id(n),iQ,iB,am)
               else
                 idd=id(n)
                 if(ist(n).eq.9)then
                   idd=sign(abs(idd)/100,idd)
                 endif
                 call idmass(idd,am)
                 if(ist(n).eq.31)
     .           am=sqrt( e**2 - px(n)**2 - py(n)**2 - pz(n)**2 )
               endif
             else
               am  = eam(n)
               e = sqrt ( px(n)**2 + py(n)**2 + pz(n)**2 + am**2 ) 
             endif
             pptl(4,nxx) = e
             pptl(5,nxx) = am
             en(n)  = e
             p5(n)  = am
             xorptl(1,nxx) = x(n)
             xorptl(2,nxx) = y(n)
             xorptl(3,nxx) = z(n)
             xorptl(4,nxx) = t(n)
             !print*,'CevtCptl ',iorptl(nxx),jorptl(nxx),'   ',nxx,'   '
             !.      ,idptl(nxx),istptl(nxx),ityptl(nxx),'   '
             !.      ,sqrt(pptl(1,nxx)**2+pptl(2,nxx)**2)
             if(n.gt.mxptl)then
               call dumpcccptl(nxx,n)
             endif
           enddo
           !------------------------------------------------------
           !write(ifch,*)'==============TreeRead==========='
           !do n=1,np
           !if(ist(n).eq.0)write(ifch,*)n,id(n),px(n),py(n),pz(n),en(n)
           !enddo
           !------------------------------------------------------
           if(ish.ge.2)call alist('list from root file&',1,np)
           call xana
           bm=bimevt
#if !__BS__ 
           call transfer_event(izmode,bm,np
     .     ,nev,npt,ngl,kol,phi,phir
     .     ,psi2,psi3,psi4,psi5,ecci2,ecci3,ecci4,ecci5
     .      ,id,ist,ity,ior,jor,px,py,pz,en,p5,x,y,z,t)
           call transfer_event_4(izmode,bm,np
     .      ,id,ist,ity,ior,jor,px,py,pz,en,p5,x,y,z,t)
#endif
           if(mod(nrevt,modsho).eq.0)
     .     write(ifmt,'(a,i7,a,a,a)')'event',
     .     nrevt,' from ',fnamein(1:index(fnamein,' ')-1),'  '
           if(hepmc.eq.1 .or. hepmc.eq.2)then
             if(ish.ge.2)call alist('Fill HepMC from Root&',0,0)
             dy=-hepmc_rapcms
c            call fillhepmc(iextree,nrevt,engy,dy,laproj,maproj
c    .       ,latarg,matarg,nhard,ngl,npartp,npartt,nspecp,nspecn,np,bm
c    .       ,sigt,id,ist,ity,ior,px,py,pz,eam,x,y,z,t
c    .       ,hepmc_record_mode,hepmc_tau_decay,hepmc_record_id_nb
c    .       ,hepmc_record_id_list)
             call fillhepmc(iextree,ihepmc3,nrevt,fegyevt,dy,laproj
     .       , maproj, latarg, matarg, nhard, ngl, npartproj, nparttarg
     .       , nspecprojp, nspecprojn, nspectargp, nspectargn, np, bm
     .       , sigt, id, ist, ity, ior, px, py, pz, eam, x, y, z, t
     .       , hepmc_record_mode, hepmc_tau_decay, hepmc_record_id_nb
     .       , hepmc_record_id_list)
           endif
         enddo
         call etreeclose()
        endif
      enddo
      if(hepmc.eq.1 .or. hepmc.eq.2) call closehepmc()
      nevent=nrevt
      call wrxx
      call swopen
      write(ifmt,'(a)')'second read in getCevtCptl'
      call aread
#if !__BS__ 
      call finish_program()
      call finish_program_4()
#endif

      write(ifmt,'(67a)')('#',i=1,24),' end getCevtCptl ',('#',i=1,26)

      end

c###################################################################################
      subroutine orderTree(maxpom)
#include "aaa.h"
      character*400 fnamemtr,fnamein, cbasout
      logical rootfile(10000)
      common/croot/imtr,fnamemtr /croot2/fnamein,ix1 
      common/croot6/cbasout,iba,nout
      common/croot4/rootfile
      common /cnnnhis/nnnhis
      parameter (mxptlxx=1.75*mxptl)
      integer id (mxptlxx), ist (mxptlxx), ity (mxptlxx)
      integer ior (mxptlxx), jor (mxptlxx)
      real px (mxptlxx), py (mxptlxx), pz (mxptlxx), en (mxptlxx)
      real p5 (mxptlxx)
      real x (mxptlxx), y (mxptlxx), z (mxptlxx), t (mxptlxx)
      common/cigrpac/igrpac /cigrval/igrval 
      common/cmxpom/mxpom
      character*200 fn
      parameter(mxifiles=200,mxpomfu=10000)
      integer npomfu(mxifiles,0:mxpomfu)
      
      mxpom=maxpom
      
      if(mxpom.gt.maxpom)stop'\n\n 21062013 \n\n'

      write(ifmt,'(67a)')('#',i=1,26),' orderTree ',('#',i=1,28)

      call getNames(2)
      call checkRootFiles(1,i1,i3,i4)
      if(i4.eq.-1)call checkRootFiles(2,i1,i3,i4)
      idel=igrTree/igrpac
      i3x=i3
      i4x=i4
      
      do kk=1,igrpac
      igrval=kk
      i3=max(i3x,i1+(kk-1)*idel)
      i4=min(i4x,i1-1+kk*idel)

      do numpom=1,mxxpom
        nfillt(numpom)=0
        nopent(numpom)=0
      !55555555555555enddo
      
      call clop(3)

      write(ifmt,'(2(a,i8),a,i3,a)')'orderTree: loop',i3,'  to',i4
     .,'                   for numpom = ', numpom

      do loo=1,2

      if(loo.eq.1)then
        !write(ifmt,'(a)')'####### First loop: count #######'
        do npom=1,mxxpom
          nfullt(npom)=0
        enddo
      else
        do  i=i3,i4
        ii3=i-i3+1
        nfull=npomfu(ii3,0)
        do nfu=1,nfull
          npom=npomfu(ii3,nfu)
          if(npom.ge.1.and.npom.le.maxpom)then
            nfullt(npom)=nfullt(npom)+1
          endif
        enddo
        enddo
        !write(ifmt,'(a)')'####### Nfullt #######'
        !do npom=1,mxxpom
        !  write(ifmt,*)npom,nfullt(npom)
        !enddo
        !write(ifmt,'(a)')'####### Second loop: store #######'
      endif

      do i=i3,i4

        ii3=i-i3+1
        call updateName(i)
        if(rootfile(i-i1+1))then
         !print*,fnamein(1:index(fnamein,' ')-1),maxpom
         call etreeopen(fnamein(1:index(fnamein,' ')-1)//CHAR(0))
         call etreenevent(nevents)
         call etreehead(iversn,laproj,maproj,latarg,matarg,engy
     .   ,nfreeze,nfull)
         if(nevents.ne.nfreeze*nfull)stop'\n\n ERROR 19062013 \n\n'
         npomfu(ii3,0)=nfull
         ii=0
         do nfu=1,nfull
         do nfr=1,nfreeze
           ii=ii+1
c needs update, not compatible with definition of etreeevt
      stop'ERROR 26092020'
c           call etreeevt(ii,np,cnt,iextree,id,ist,ity,ior,jor
c     .     ,px,py,pz,en,x,y,z,t)
           nptl=np
           npom=0
           do n=1,np
           call checkcccptl(n)
           nxx=n
           if(n.gt.mxptl)nxx=mxptl+1
           idptl(nxx)    = id(n)
           istptl(nxx)   = ist(n)
           ityptl(nxx)   = ity(n)
           iorptl(nxx)   = ior(n)
           jorptl(nxx)   = jor(n)
           pptl(1,nxx)   = px(n)
           pptl(2,nxx)   = py(n)
           pptl(3,nxx)   = pz(n)
           pptl(4,nxx)   = en(n)
           if(ity(n).eq.61)then
             call getQBM(id(n),iQ,iB,am)
           else
             call idmass(id(n),am)
           endif
           pptl(5,nxx)   = am
           p5(n)  = am
           xorptl(1,nxx) = x(n)
           xorptl(2,nxx) = y(n)
           xorptl(3,nxx) = z(n)
           xorptl(4,nxx) = t(n)
           if(n.gt.mxptl)call dumpcccptl(nxx,n)
           if(ist(n).eq.30.or.ist(n).eq.31)npom=npom+1
           enddo
           if(loo.eq.1)then
             if(izmode.eq.6)then !number of Pomerons 
               continue !  default
             elseif(izmode.eq.7)then !vzero
               npom=centrVar(7)
             else
               stop'\n\n ERROR 23062013 \n\n'
             endif
             if(nfull.gt.mxpomfu)stop'\n\n 24062013\n\n'
             if(ii3.gt.mxifiles)stop'\n\n 24062013b\n\n'
             if(nfr.eq.1)then
               npomfu(ii3,nfu)=npom
             else
               npomfu(ii3,nfu)=npomfu(ii3,nfu)+npom
             endif
             if(nfr.eq.nfreeze)npomfu(ii3,nfu)
     .                  =nint(npomfu(ii3,nfu)/float(nfreeze))
           else
             npom=npomfu(ii3,nfu)
           endif
           if(loo.eq.2)then
            if(npom.eq.numpom)then !55555555555555555
             if(npom.ge.1.and.npom.le.maxpom)then
               ! write(ifmt,*)fnamein(1:index(fnamein,' ')-1)
               !.,npom,nfullt(npom)
               if(nfullt(npom).ge.2)
     .         call treestore(npom)
             endif
            endif!55555555555555555555555
           endif 
         enddo
         enddo
         call etreeclose()
        endif
        
      enddo ! i

      enddo  ! loo
      
      enddo ! numpom 55555555555555

      write(ifmt,'(a,$)')'closetree'
      do npom=1,mxxpom
       if(nopent(npom).gt.0)then
        write(ifmt,'(i3,$)')npom
        call tfname(npom,fn,nfn)          
        call closetree(npom,iextree,ihepmc3,fn(1:nfn)//CHAR(0))
       endif
      enddo
      write(ifmt,'(a)')' '
       
      call clop(3)

      enddo ! kk

      write(ifmt,'(67a)')('#',i=1,24),' end orderTree ',('#',i=1,26)
      stop'end of orderTree'

      end

c###################################################################################
      function fmux(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
      !   xpara(1,n) ... etamin
      !   xpara(2,n) ... etamax
      !   xpara(3,n) ... ptmin
      !   xpara(4,n) ... ptmax
      !   xpara(5,n) ... factor
      !   xpara(6,n) ... divisor
      !   xpara(7,n) ... absolute value of eta (1)
      !   xpara(8,n) ... yboost if >< 0.0001 (boost back)
      !   xpara(9,n) ... etamin2 if >< 0.
      !   xpara(10,n) ... etamax2 if >< 0.
      !   xpara(11,n) ... x0 = normal
      !                   x2 = percentile  
      !                   0x = multiplicity
      !                   2x = Et
#include "xan.h"
      n=1
      xpara(1,n)=p1
      xpara(2,n)=p2
      xpara(3,n)=p3
      xpara(4,n)=p4
      xpara(5,n)=p5
      xpara(6,n)=p6
      xpara(7,n)=p7
      xpara(8,n)=p8
      xpara(9,n)=p9
      xpara(10,n)=p10
      xpara(11,n)=p11
      call mux(n,0)
      fmux=ypara(1,n)
      end

c###################################################################################
      subroutine getFsiParameters(i1,i2,i4,i5,i6)
      integer*4 FsiParameters(15)
      common/cFsiParameters/FsiParameters

      FsiParameters(1)=i1  !particle 1
      FsiParameters(2)=i2  !particle 2
      FsiParameters(4)=i4  !switch coulomb interaction on(1)/off(0)
      FsiParameters(5)=i5  !switch quantum statistics on(1)/off(0)
      FsiParameters(6)=i6  !switch strong interaction on(1)/off(0)

      FsiParameters(3)=1    !if set to 0 default fsi parameters will be used
      FsiParameters(7)=0    !switch 3rd body Coulomb influence on(1)/off(0)
      FsiParameters(8)=2    !maproj 1 + matarg 1 ->mass of the third body
      FsiParameters(9)=2    !laproj 1 + latarg 1 ->charge of the third body
      FsiParameters(10)=1   !sign of the 3rd body charge(+1); -1 for test
      FsiParameters(11)=1   !1 use spherical approximation
      FsiParameters(12)=0   !0 use square well approximation
      FsiParameters(13)=0   !--not used yet
      FsiParameters(14)=0   !--not used yet
      FsiParameters(15)=0   !fastQSon(1)/off(0): 0=Lednicky FSI weight
                            !1=else cos(delta_p*delta_x)
      end



c#######################################################################################
c#######################################################################################
c############                       NOT BASIC
c#######################################################################################
c#######################################################################################


#if !__BS__ 


      subroutine d2hstore(nev,nfr)
#include "aaa.h"
      common/cen/ncentr
      common/cranphi/ranphi
      double precision etamax
      parameter (mxptlyy=mxptl/20,nxbins=40,nybins=20, etamax=2d0)
      parameter (mxfr=10)
      real eta(2,mxfr,mxptlyy),phi(2,mxfr,mxptlyy),pta(2,mxfr,mxptlyy)
      integer nmb(2,mxfr)
      double precision deleta, delphi
      parameter                          (kmax=3)
      logical go,goi,goj,gok(kmax)
      data nfilld2h /0/
      save nfilld2h,nmb
      if(nfr+1.gt.mxfr)stop'\n\n STOP in d2hstore, mxfr too small \n\n'

      do k=1,kmax
      gok(k)=.false.
      enddo
      do k=1,kmax
      b1=rclass(2,1,1)   !temporalily
      b2=rclass(2,2,1)   !temporalily
      if(bimevt.ge.b1.and.bimevt.le.b2)gok(k)=.true.
      enddo
      go=.false.
      do k=1,kmax
      if(gok(k))go=.true.
      enddo
      if(.not.go)return

      ne=1+mod(nev-1,2)
      n=0
      do i=nptlpt+1,nptl
        if(istptl(i).eq.0)then
          pt=sqrt(  pptl(1,i)**2 + pptl(2,i)**2 )
          pz=pptl(3,i)
          if(pz.ne.0..and.pt.ge.2.)then
          n=n+1
          pta(ne,nfr+1,n)=pt
          eta(ne,nfr+1,n)=sign(1.,pz)*
     *       alog((sqrt(pz**2+pt**2)+abs(pz))/pt)
          p1=pptl(1,i)
          p2=pptl(2,i)
          !~~~~~
          phinll=phievt+ranphi
          aa=cos(phinll)
          bb=sin(phinll)
          cc=-sin(phinll)
          dd=cos(phinll)
          px=p1*aa+p2*bb
          py=p1*cc+p2*dd
          !~~~~~~
          phi(ne,nfr+1,n)=sign(1.,py)*acos(px/pt)
          endif
        endif
      enddo
      nmb(ne,nfr+1)=n

      if(mod(nev,2).ne.0)return
      if(mod(nfr+1,nfreeze).ne.0)return

      nfilld2h=nfilld2h+1
      if(nfilld2h.eq.1)then
        nopend2h=nopend2h+1
        ii=index(fndt(1:nfndt),".data")-1
        write(ifmt,'(a)')'open2dhisto '//fndt(1:ii)//'.d2h'
        call open2dhisto(fndt(1:ii)//'.d2h'//CHAR(0)
     .  ,kmax,nxbins,nybins,etamax)
      endif

      !---same event---
      do kk=1,kmax
      if(gok(kk))then
      do nn=1,2
      do nf=1,nfreeze
      n=nmb(nn,nf)
      do i=1,n
      pti=pta(nn,nf,i)
      if(kk.eq.1)goi=   pti.ge.3. .and. pti.le.4.
      if(kk.eq.2)goi=   pti.ge.4. .and. pti.le.6.
      if(kk.eq.3)goi=   pti.ge.2. .and. pti.le.3.
      if(goi)then
      do j=1,n
      ptj=pta(nn,nf,j)
      if(kk.eq.1)goj=   ptj.ge.2. .and. ptj.le.pti
      if(kk.eq.2)goj=   ptj.ge.2. .and. ptj.le.pti
      if(kk.eq.3)goj=   ptj.ge.2. .and. ptj.le.3.
      if(j.ne.i.and.goj)then
      deleta=eta(nn,nf,i)-eta(nn,nf,j)
      delphi=phi(nn,nf,i)-phi(nn,nf,j)
      if(delphi.gt.pi)delphi=delphi-2*pi
      if(delphi.lt.-pi)delphi=delphi+2*pi
      call fill2dhisto1(kk,deleta,delphi)
      endif
      enddo
      endif
      enddo
      enddo
      enddo
      endif
      enddo

      !---mixed---
      do kk=1,3
      if(gok(kk))then
      do nf1=1,nfreeze
      n1=nmb(1,nf1)
      do nf2=1,nfreeze
      n2=nmb(2,nf2)
      do i=1,n1
      pti=pta(1,nf1,i)
      if(kk.eq.1)goi=   pti.ge.3. .and. pti.le.4.
      if(kk.eq.2)goi=   pti.ge.4. .and. pti.le.6.
      if(kk.eq.3)goi=   pti.ge.2. .and. pti.le.3.
      if(goi)then
      do j=1,n2
      ptj=pta(2,nf2,j)
      if(kk.eq.1)goj=   ptj.ge.2. .and. ptj.le.pti
      if(kk.eq.2)goj=   ptj.ge.2. .and. ptj.le.pti
      if(kk.eq.3)goj=   ptj.ge.2. .and. ptj.le.3.
      if(goj)then
      deleta=eta(1,nf1,i)-eta(2,nf2,j)
      delphi=phi(1,nf1,i)-phi(2,nf2,j)
      if(delphi.gt.pi)delphi=delphi-2*pi
      if(delphi.lt.-pi)delphi=delphi+2*pi
      call fill2dhisto2(kk,deleta,delphi)
      endif
      enddo
      endif
      enddo
      enddo
      enddo
      endif
      enddo

      end


      subroutine d2hclose
#include "aaa.h"
      if(nopend2h.gt.0)then
        write(ifmt,'(a)')'close2dhisto'
        call close2dhisto()
      endif
      end

cccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine dih   !removed in 3444m
c   also removed: dih.cpp dih.h 
c    in KWb/ and in Makefile
cccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccc
c  subroutine femto !***** removed in 3444m4
c   also removed: fto.cpp fto.h fsiwled.f fsitools.f
c    in KWb/ and in Makefile
cccccccccccccccccccccccccccccc

      subroutine vmd
#include "aaa.h"
      real cMin(10),   cMax(10)
      character*400 fnamemtr
      common/croot/imtr,fnamemtr
      common/cNpart/sumNpart
      k1=ivmd
      call getNames(2)
      iMaxNpart=sumNpart/nfull/nfreeze
      nClass=nrclass(k1)   !number of centrality classes (<=10)
      do i=1,nrclass(k1)
      cMin(i)  =   rclass(k1,1,i)
      cMax(i)  =   rclass(k1,2,i)
      !print*,i,cMin(i),cMax(i)
      enddo

      write(ifmt,'(a)')'enter eposvmd'
      call clop(2)

c      call eposvmdm(
c     . fnamein(1:index(fnamein,' ')-1)//CHAR(0),
c     . cbasout(1:index(cbasout,' ')-1)//CHAR(0),
c     . nout,
c     . fnamemtr(1:index(fnamemtr,' ')-1)//CHAR(0),
c     . nClass, cMin, cMax,
c     . mixevt, iMaxNpart
c     . )

      call clop(1)
      write(ifmt,'(a)')'exit eposvmd'

      end


#endif

