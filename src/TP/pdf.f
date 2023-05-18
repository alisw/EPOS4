C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

csp------------------------------------------------------------------------
csp     This file contains everything related to PDF
csp------------------------------------------------------------------------

ctp------------------------------------------------------------------
      function psdpdf(xxx,qqs,qq,icq,iq)
ctp------------------------------------------------------------------
ctp PDF used every where in EPOS
ctp iq>0 = val quark (if exist)
ctp iq<0 = sea quark
ctp icq : hadron type
ctp------------------------------------------------------------------
#include "sem.h"

      x=xxx!  min(xxx,0.99999)     !needed for idrafl

      kref=icq
      psdpdf=psdfh4(x,qqs,qq,icq,iq)   !GRV94

c alternative CTEQ 6 (new)
c      kref=2          !only proton
c      psdpdf=cteq6(x,qqs,qq,icq,iq)

c alternative GRV 2008 (new)
c      kref=2          !only proton
c      psdpdf=grv08(x,qqs,qq,icq,iq) !problems at high energy


      return
      end

c------------------------------------------------------------------------
      function psjpdf(pdfunc,qt,s,y0)
c-----------------------------------------------------------------------
c Integral over y of pdf = "pdfunc"
c (former  psjgrv081, psjcteq61, psjvrg1 and psjwo1) and correspond to
c ffsigi (cross-section for the hard processes (dijet))
c-----------------------------------------------------------------------
      common /ar3/   x1(7),a1(7)
      common /cnsta/ pi,pii,hquer,prom,piom,ainfin
#include "sem.h"
      double precision xt,ymin,ymax,y,xmin,xmax,xx1,xx2,sh,z,ft,fx,t
     *,pdfunc,psjpdfd
      external pdfunc


      psjpdf=0.
      if(s.le.4.*qt)return
      psjpdfd=0.d0

      xt=2.d0*sqrt(dble(qt)/dble(s))
      ymax=min(dble(y0),log(1d0/xt+sqrt((1d0/xt-1d0)*(1d0/xt+1d0))))
      ymin=-ymax

      do i=1,7
      do m=1,2
        y=.5d0*(ymax+ymin+(ymin-ymax)*dble((2*m-3)*x1(i)))
        xmin=xt**2/2.d0/(2.d0-xt*exp(-y))
        xmax=1.d0-xt*exp(y)/2.d0

        fx=0.d0
        do i1=1,7
        do m1=1,2
          xx1=xt*exp(y)/2d0+xmin*(xmax/xmin)**dble(.5+x1(i1)*(m1-1.5))
          xx2=xt*exp(-y)*xx1/(2.d0*xx1-xt*exp(y))
          z=xx1*xx2
          sh=z*s
          t=sh/2d0*(1d0-sqrt(max(0d0,1d0-4d0*dble(qt)/sh)))
          ft=pdfunc(t,qt,sngl(xx1),sngl(xx2),sh)
          fx=fx+dble(a1(i1))*ft/sh**2
        enddo
        enddo
        fx=fx*0.5d0*log(xmax/xmin)
        psjpdfd= psjpdfd+dble(a1(i))*fx
      enddo
      enddo
       psjpdf= sngl(psjpdfd*(ymax-ymin))*pi**3
     **pssalf(qt/qcdlam)**2*sqrt(qt)
     **2                            !2 jets are produced per hard process
      return
      end

csp**************************************************************
csp                    GJR 08 (GRV2008)
csp**************************************************************

c-----------------------------------------------------------------------
      double precision function psjgrv08x(t,qt,xx1,xx2,s)
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s,g1,ub1,u1,db1,d1,sb1,s1
     *                           ,g2,ub2,u2,db2,d2,sb2,s2

      g1=dble(grv08(xx1,qt,0.,2,0))
      ub1=dble(grv08(xx1,qt,0.,2,-1))
      u1=dble(grv08(xx1,qt,0.,2,1))+ub1
      db1=dble(grv08(xx1,qt,0.,2,-2))
      d1=dble(grv08(xx1,qt,0.,2,2))+db1
      sb1=dble(grv08(xx1,qt,0.,2,-3))
      s1=sb1
      g2=dble(grv08(xx2,qt,0.,2,0))
      ub2=dble(grv08(xx2,qt,0.,2,-1))
      u2=dble(grv08(xx2,qt,0.,2,1))+ub2
      db2=dble(grv08(xx2,qt,0.,2,-2))
      d2=dble(grv08(xx2,qt,0.,2,2))+db2
      sb2=dble(grv08(xx2,qt,0.,2,-3))
      s2=sb2

      psjgrv08x=g1*g2*(psbori(99,s,t,0,0,1)+psbori(99,s,s-t,0,0,1)
     *+psbori(99,s,t,0,0,2)+psbori(99,s,s-t,0,0,2))/2. !bg

     *+(psbori(99,s,t,0,1,1)+psbori(99,s,s-t,0,1,1))*
     *(g2*(u1+ub1+d1+db1+s1+sb1)+g1*(u2+ub2+d2+db2+s2+sb2))

     *+(psbori(99,s,t,1,1,1)+psbori(99,s,s-t,1,1,1))/2.* !bg
     *(u1*u2+ub1*ub2+d1*d2+db1*db2+s1*s2+sb1*sb2)

     *+(psbori(99,s,t,1,-1,1)+psbori(99,s,s-t,1,-1,1)
     *+psbori(99,s,t,1,-1,2)+
     *psbori(99,s,s-t,1,-1,2)+psbori(99,s,t,1,-1,3)
     *+psbori(99,s,s-t,1,-1,3))*
     *(u1*ub2+ub1*u2+d1*db2+db1*d2+s1*sb2+sb1*s2)

     *+(psbori(99,s,t,1,2,1)+psbori(99,s,s-t,1,2,1))*
     *((u1+ub1)*(d2+db2+s2+sb2)+(u2+ub2)*(d1+db1+s1+sb1)+
     *(d1+db1)*(u2+ub2+s2+sb2)+(d2+db2)*(u1+ub1+s1+sb1)+
     *(s1+sb1)*(u2+ub2+d2+db2)+(s2+sb2)*(u1+ub1+d1+db1))

      return
      end



csp----------------------------------------------------------------------
          function grv08(x,qqs,qq,icq,iq)
csp----------------------------------------------------------------------
csp GRV 2008 pdf
csp----------------------------------------------------------------------
      double precision xuv(118,99),xdv(118,99),xgl(118,99),xub(118,99),
     *xdb(118,99),xsb(118,99),
     *grid(217),arg(2),dfint
      integer ng(2)
       common /gridc/ grid,ng
       common /xuvc/ xuv
        common /xdvc/ xdv
        common /xglc/ xgl
        common /xubc/ xub
        common /xdbc/ xdb
        common /xsbc/ xsb
         common/cpigrv08/npigrv08
        data npigrv08/0/
        npigrv08=npigrv08+1
      if(npigrv08.eq.1)call GJR08FFNSinit(15)


       dum=qq
       idum=icq

      if(x.gt..99999)then
        grv08=0.
        return
      endif
      arg(1)=x
      arg(2)=qqs


       if(iq.eq.0)then
       grv08=dfint(2,arg,ng,grid,xgl)
       elseif(iq.eq.1)then
       grv08=dfint(2,arg,ng,grid,xuv)
       elseif(iq.eq.2)then
       grv08=dfint(2,arg,ng,grid,xdv)
       elseif(iq.eq.-3)then
       grv08=dfint(2,arg,ng,grid,xsb)
       elseif(iq.eq.-1)then
       grv08=dfint(2,arg,ng,grid,xub)
       elseif(iq.eq.-2)then
       grv08=dfint(2,arg,ng,grid,xdb)
       else
       grv08=0.
       endif

        end


csp----------------------------------------------------------------------
!! CERNLIB E104 modified to be used with GJR08 GRIDS:
!! Name changed from fint to dfint.
!! Real variables changed to double precision.
!! External references to CERNLIB (error handling) routines removed.
          DOUBLE PRECISION FUNCTION DFINT(NARG,ARG,NENT,ENT,TABLE)
          INTEGER   NENT(*), INDEX(32)
          DOUBLE PRECISION ARG(*),   ENT(*),TABLE(*), WEIGHT(300)
          DFINT  =  0d0
          IF(NARG .LT. 1  .OR.  NARG .GT. 5)  GOTO 300
          LMAX      =  0
          ISTEP     =  1
          KNOTS     =  1
          INDEX(1)  =  1
          WEIGHT(1) =  1d0
          DO 100    N  =  1, NARG
             X     =  ARG(N)
             NDIM  =  NENT(N)
             LOCA  =  LMAX
             LMIN  =  LMAX + 1
             LMAX  =  LMAX + NDIM
             IF(NDIM .GT. 2)  GOTO 10
             IF(NDIM .EQ. 1)  GOTO 100
             H  =  X - ENT(LMIN)
             IF(H .EQ. 0.)  GOTO 90
             ISHIFT  =  ISTEP
             IF(X-ENT(LMIN+1) .EQ. 0d0)  GOTO 21
             ISHIFT  =  0
             ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
             GOTO 30
  10         LOCB  =  LMAX + 1
  11         LOCC  =  (LOCA+LOCB) / 2
             IF(X-ENT(LOCC))  12, 20, 13
  12         LOCB  =  LOCC
             GOTO 14
  13         LOCA  =  LOCC
  14         IF(LOCB-LOCA .GT. 1)  GOTO 11
             LOCA    =  MIN0( MAX0(LOCA,LMIN), LMAX-1 )
             ISHIFT  =  (LOCA - LMIN) * ISTEP
             ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
             GOTO 30
  20         ISHIFT  =  (LOCC - LMIN) * ISTEP
  21         DO 22  K  =  1, KNOTS
                INDEX(K)  =  INDEX(K) + ISHIFT
  22            CONTINUE
             GOTO 90
  30         DO 31  K  =  1, KNOTS
                INDEX(K)         =  INDEX(K) + ISHIFT
                INDEX(K+KNOTS)   =  INDEX(K) + ISTEP
                WEIGHT(K+KNOTS)  =  WEIGHT(K) * ETA
                WEIGHT(K)        =  WEIGHT(K) - WEIGHT(K+KNOTS)
  31            CONTINUE
             KNOTS  =  2*KNOTS
  90         ISTEP  =  ISTEP * NDIM
 100         CONTINUE
          DO 200    K  =  1, KNOTS
             I  =  INDEX(K)
             DFINT  =  DFINT + WEIGHT(K) * TABLE(I)
 200         CONTINUE
          RETURN
 300      WRITE(*,1000) NARG
          STOP
1000      FORMAT( 7X, 24HFUNCTION DFINT... NARG =,I6,
     +              17H NOT WITHIN RANGE)
          END

csp----------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GJR08FFNS GRIDS (Eur. Phys. J. C53 (2008) 355 and arXiv:0709.0614) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This package contains the GJR08FFNSFFNS LO, NLO(MSbar) and NLO(DIS)
!! dynamical parton distributions of the nucleon and their associated
!! exact iterative solutions for alpha_s. The heavy quark masses used as
!! well as the lambda_QCD parameters for the approximate 'asymptotic'
!! solution for alpha_s are given in the paper.
!!
!! The sets resulting from displacements in the parameter space with
!! respect to the central value of the NLO(MSbar) fit along the plus
!! (minus) directions of the (rescaled) eigenvectors of the hessian
!! matrix are included. The tolerance parameter for these displacements
!! was chosen to be T=4.654 for a total of 1739 fitted points. Since
!! alpha_s was included as a free parameter in the error estimation the
!! use of the provided alpha_s solution for each set is mandatory for
!! uncertainty studies (the difference on alpha_s for different
!! eigenvector sets can be as large as 10% at low scales).
!!
!! The predictions for F2_light (including target mass corrections),
!! F2_charm and F2_bottom DIS structure functions of the proton arising
!! from our NLO(MSbar) dynamical distributions are also included.
!!
!! The grids are generated for the following ranges:
!!  parton densities:    10^-9 <= x <= 1 and Qo^2 <= Q^2 <= 10^8 (GeV^2)
!!  alpha_s:                                 Qo^2 <= Q^2 <= 10^8 (GeV^2)
!!  structure functions: 10^-9 <= x <= 1 and Qo^2 <= Q^2 <= 10^5 (GeV^2)
!! where Qo^2 = 0.5 GeV^2 (0.3 GeV^2) for the NLO (LO) fits. Outside
!! these ranges the output is either set to NaN or obtained through
!! extrapolation and should NOT be used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The package has been written in standard Fortran77 and tested using
!! gfortran(v4.2.1) and ifort(v9.1.043) and g77(v3.4.6). The routines
!! use a modification of the standard multidimensional linear
!! interpolation routine FINT (CERNLIB E104) distributed as the file
!! 'dfint.f'. A mimimum of about 18.4 MB of memory is needed to load the
!! grid. For questions, comments, problems etc please contact:
!! pjimenez@het.physik.uni-dortmund.de
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GJR08FFNSinit(set):
!!      Initialization routine of the package to be called before using
!!      any of the other subroutines. It reads the file
!!      './GJR08FFNS.grd', where ./ may be replaced by the path from the
!!      working directory to the grid file.
!!        set == Index specifying the set to be loaded:
!!             0 NLO(MSbar).
!!             1,2,...,13 set corresponding to a displacement +T with
!!                        respect to the set 0 in the direction of the
!!                        ith eigenvector.
!!             -1,-2,...,-13 the same for displacements -T.
!!             14 NLO(DIS)
!!             15 LO
!! GJR08FFNSx'parton'(x,Q2) with 'parton' = uv,dv,gl,ub,db,sb:
!!      Parton distribution 'parton' (times x) for the set specified in
!!      the last call to GJR08FFNSinit.
!!        x == Bjorken-x.
!!        Q2 == Q**2 (GeV**2) == Factorization scale
!!                            == Renormalization scale.
!! GJR08FFNSalphas(Q2):
!!      Value of alpha_s (no additional 2pi or 4pi factors.) for the set
!!      specified in the last call to GJR08FFNSinit.
!!        Q2 == Q**2 (GeV**2) == Renormalization scale.
!! GJR08FFNS'function'(x,Q2)
!!                with 'function' = F2l(light), F2c(charm), F2b(bottom):
!!      Proton structure function 'function' for the NLO(MSbar) fit if
!!      set was 0 in the last call to GJR08FFNSinit, NaN otherwise.
!!        x == Bjorken-x
!!        Q2 == Q**2 (GeV**2) == DIS Momentum transfer squared.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GJR08FFNSinit(set)
       implicit none
      integer mxho
      parameter(mxho=10)
      character*500 fnch,fnhi,fndt,fnhm,fnii,fnid,fnie,fnrj,fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnio,fnho,fn3g,fn3p,fn3d
     *,fn3f,fn3f1,fn3f2,fn3f3,fn3f4,fn3f5,fnhpf
      common/fname/ fnch,fnhi,fndt,fnhm,fnii,fnid,fnie,fnrj,fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnio,fnho(mxho),fn3g,fn3p,fn3d
     *,fn3f,fn3f1,fn3f2,fn3f3,fn3f4,fn3f5,fnhpf
      integer nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnio,nfnho,nfn3g,nfnhm
     *,nfn3p,nfn3d,nfn3f,nfn3f1,nfn3f2,nfn3f3,nfn3f4,nfn3f5,nfnhpf
      common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnio,nfnho(mxho),nfn3g,nfnhm
     *,nfn3p,nfn3d,nfn3f,nfn3f1,nfn3f2,nfn3f3,nfn3f4,nfn3f5,nfnhpf
       integer ng(2),init,set,i,j,k,l
       double precision fgrid(-13:16,-3:3,118,99),grid(217),
     &                   xuv(118,99),xdv(118,99),xgl(118,99),
     &                   xub(118,99),xdb(118,99),xsb(118,99),
     &                   alphas(118,99),F2l(118,99),F2c(118,99),
     &                   F2b(118,99)
       save init,fgrid
       common /gridc/ grid,ng
       common /xuvc/ xuv
       common /xdvc/ xdv
       common /xglc/ xgl
       common /xubc/ xub
       common /xdbc/ xdb
       common /xsbc/ xsb
       common /alphasc/ alphas
       common /F2lc/ F2l
       common /F2cc/ F2c
       common /F2bc/ F2b
       if (init/=1) then
        init=1
        data ng /118,99/
        data grid
     &   /1d-9,1.25d-9,1.6d-9,2d-9,2.5d-9,3.16d-9,4d-9,5d-9,6.3d-9,8d-9,
     &    1d-8,1.25d-8,1.6d-8,2d-8,2.5d-8,3.16d-8,4d-8,5d-8,6.3d-8,8d-8,
     &    1d-7,1.25d-7,1.6d-7,2d-7,2.5d-7,3.16d-7,4d-7,5d-7,6.3d-7,8d-7,
     &    1d-6,1.25d-6,1.6d-6,2d-6,2.5d-6,3.16d-6,4d-6,5d-6,6.3d-6,8d-6,
     &    1d-5,1.25d-5,1.6d-5,2d-5,2.5d-5,3.16d-5,4d-5,5d-5,6.3d-5,8d-5,
     &    1d-4,1.25d-4,1.6d-4,2d-4,2.5d-4,3.16d-4,4d-4,5d-4,6.3d-4,8d-4,
     &    1d-3,1.25d-3,1.6d-3,2d-3,2.5d-3,3.16d-3,4d-3,5d-3,6.3d-3,8d-3,
     &    1d-2,1.25d-2,1.6d-2,2d-2,2.5d-2,3.16d-2,4d-2,5d-2,6.3d-2,8d-2,
     &    0.10d0,0.125d0,0.15d0,0.175d0,0.20d0,0.225d0,0.25d0,0.275d0,
     &    0.30d0,0.325d0,0.35d0,0.375d0,0.40d0,0.425d0,0.45d0,0.475d0,
     &    0.50d0,0.525d0,0.55d0,0.575d0,0.60d0,0.625d0,0.65d0,0.675d0,
     &    0.70d0,0.725d0,0.75d0,0.775d0,0.80d0,0.825d0,0.85d0,0.875d0,
     &    0.9d0,0.920d0,0.94d0,0.960d0,0.98d0,1d0,
     &    0.3d0,0.31d0,0.35d0,0.375d0,0.4d0,0.45d0,0.5d0,0.51d0,0.525d0,
     &    0.55d0,0.575d0,0.6d0,0.65d0,0.7d0,0.75d0,0.8d0,0.85d0,0.9d0,
     &    1d0,1.25d0,1.6d0,2d0,2.5d0,3.16d0,4d0,5d0,6.3d0,8d0,
     &    1d1,1.25d1,1.6d1,2d1,2.5d1,3.16d1,4d1,5d1,6.3d1,8d1,
     &    1d2,1.25d2,1.6d2,2d2,2.5d2,3.16d2,4d2,5d2,6.3d2,8d2,
     &    1d3,1.25d3,1.6d3,2d3,2.5d3,3.16d3,4d3,5d3,6.3d3,8d3,
     &    1d4,1.25d4,1.6d4,2d4,2.5d4,3.16d4,4d4,5d4,6.3d4,8d4,
     &    1d5,1.25d5,1.6d5,2d5,2.5d5,3.16d5,4d5,5d5,6.3d5,8d5,
     &    1d6,1.25d6,1.6d6,2d6,2.5d6,3.16d6,4d6,5d6,6.3d6,8d6,
     &    1d7,1.25d7,1.6d7,2d7,2.5d7,3.16d7,4d7,5d7,6.3d7,8d7,1d8/
        open(10,file=fngrv(1:nfngrv)//'GJR08FFNS.grd')
        do 4 i=-13,13
         do 3 j=1,118
          do 2 k=1,99
           read(10,*) fgrid(i,-3,j,k),fgrid(i,-2,j,k),fgrid(i,-1,j,k),
     &                fgrid(i,0,j,k),fgrid(i,1,j,k),fgrid(i,2,j,k),
     &                fgrid(i,3,j,k)
           do 1 l=-3,3
            if (grid(118+k) < 0.5d0) then
             fgrid(i,l,j,k)=0d0
csp             fgrid(i,l,j,k)=0d0/fgrid(i,l,j,k)
            end if
    1      continue
    2     continue
    3    continue
    4   continue
        i=14
         do 7 j=1,118
          do 6 k=1,99
           read(10,*) fgrid(i,-3,j,k),fgrid(i,-2,j,k),fgrid(i,-1,j,k),
     &                fgrid(i,0,j,k),fgrid(i,1,j,k),fgrid(i,2,j,k),
     &                fgrid(i,3,j,k)
           do 5 l=-3,3
            if (grid(118+k) < 0.5d0) then
             fgrid(i,l,j,k)=0d0
csp             fgrid(i,l,j,k)=0d0/fgrid(i,l,j,k)
            end if
    5      continue
    6     continue
    7    continue

        i=15
         do 10 j=1,118
          do 9 k=1,99
           read(10,*) fgrid(i,-3,j,k),fgrid(i,-2,j,k),fgrid(i,-1,j,k),
     &                fgrid(i,0,j,k),fgrid(i,1,j,k),fgrid(i,2,j,k),
     &                fgrid(i,3,j,k)
    9     continue
   10    continue
        i=16
         do 13 j=1,118
          do 12 k=1,99
           read(10,*) fgrid(i,1,j,k),fgrid(i,2,j,k),fgrid(i,3,j,k)
           do 11 l=1,3
            if ((grid(118+k) < 0.5d0).or.(grid(118+k) > 1d5)) then
             fgrid(i,l,j,k)=0d0
csp             fgrid(i,l,j,k)=0d0/fgrid(i,l,j,k)
            end if
   11      continue
   12     continue
   13    continue
        close(10)

       end if
       do 15 j=1,118
        do 14 k=1,99
         xuv(j,k) = fgrid(set,1,j,k)
         xdv(j,k) = fgrid(set,2,j,k)
         xgl(j,k) = fgrid(set,0,j,k)
         xub(j,k) = fgrid(set,-1,j,k)
         xdb(j,k) = fgrid(set,-2,j,k)
         xsb(j,k) = fgrid(set,-3,j,k)
         alphas(j,k) = fgrid(set,3,j,k)
         if (set==0) then
          F2l(j,k) = fgrid(16,1,j,k)
          F2c(j,k) = fgrid(16,2,j,k)
          F2b(j,k) = fgrid(16,3,j,k)
         else
          F2l(j,k) = 0d0
csp          F2l(j,k) = 0d0/F2l(j,k)
          F2c(j,k) = 0d0
csp          F2c(j,k) = 0d0/F2c(j,k)
          F2b(j,k) = 0d0
csp          F2b(j,k) = 0d0/F2b(j,k)
         end if
   14   continue
   15  continue
       return
      end  ! subroutine GJR08FFNSinit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxuv(x,Q2)
       implicit none
       integer ng(2)
       double precision xuv(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xuvc/ xuv
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxuv = dfint(2,arg,ng,grid,xuv)
       return
      end ! function GJR08FFNSxuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxdv(x,Q2)
       implicit none
       integer ng(2)
       double precision xdv(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xdvc/ xdv
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxdv = dfint(2,arg,ng,grid,xdv)
       return
      end ! function GJR08FFNSxdv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxgl(x,Q2)
       implicit none
       integer ng(2)
       double precision xgl(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xglc/ xgl
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxgl = dfint(2,arg,ng,grid,xgl)
       return
      end !function GJR08FFNSxgl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxub(x,Q2)
       implicit none
       integer ng(2)
       double precision xub(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xubc/ xub
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxub = dfint(2,arg,ng,grid,xub)
       return
      end !function GJR08FFNSxub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxdb(x,Q2)
       implicit none
       integer ng(2)
       double precision xdb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xdbc/ xdb
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxdb = dfint(2,arg,ng,grid,xdb)
       return
      end !function GJR08FFNSxdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSxsb(x,Q2)
       implicit none
       integer ng(2)
       double precision xsb(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /xsbc/ xsb
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSxsb = dfint(2,arg,ng,grid,xsb)
       return
      end !function GJR08FFNSxsb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSalphas(Q2)
       implicit none
       integer ng(2)
       double precision alphas(118,99),grid(217),arg(2),Q2,dfint
       common /gridc/ grid,ng
       common /alphasc/ alphas
       arg(1) = 1d-9
       arg(2) = Q2
       GJR08FFNSalphas = dfint(2,arg,ng,grid,alphas)
       return
      end !function GJR08FFNSalphas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSF2l(x,Q2)
       implicit none
       integer ng(2)
       double precision F2l(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /F2lc/ F2l
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSF2l = dfint(2,arg,ng,grid,F2l)
       return
      end !function GJR08FFNSF2l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSF2c(x,Q2)
       implicit none
       integer ng(2)
       double precision F2c(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /F2cc/ F2c
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSF2c = dfint(2,arg,ng,grid,F2c)
       return
      end !function GJR08FFNSF2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function GJR08FFNSF2b(x,Q2)
       implicit none
       integer ng(2)
       double precision F2b(118,99),grid(217),arg(2),x,Q2,dfint
       common /gridc/ grid,ng
       common /F2bc/ F2b
       arg(1) = x
       arg(2) = Q2
       GJR08FFNSF2b = dfint(2,arg,ng,grid,F2b)
       return
      end !function GJR08FFNSF2b

csp--------------------------------------------------------------------
csp--------------------------------------------------------------------
csp                                CTEQ6
csp--------------------------------------------------------------------
csp--------------------------------------------------------------------

C============================================================================
C                CTEQ Parton Distribution Functions: version 6.0-6.6
C                             April 10, 2002, v6.01
C                             February 23, 2003, v6.1
C                             August 6, 2003, v6.11
C                             December 12, 2004, v6.12
C                             December 4, 2006, v6.5 (CTEQ6.5M series added)
C                             March 23, 2007, v6.51 (CTEQ6.5S/C series added)
C                             April 24, 2007, v6.52 (minor improvement)
C                             March 30, 2008, v6.6
C
C   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
C       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
C       JHEP 0207:012(2002), hep-ph/0201195
C
C   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
C       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
C       JHEP 0310:046(2003), hep-ph/0303013
C
C   Ref[3]: "Neutrino dimuon Production and Strangeness Asymmetry of the Nucleon"
C       By: F. Olness, J. Pumplin, S. Stump, J. Huston, P. Nadolsky, H.L. Lai, S. Kretzer, J.F. Owens, W.K. Tung
C       Eur. Phys. J. C40:145(2005), hep-ph/0312323
C
C   Ref[4]: "CTEQ6 Parton Distributions with Heavy Quark Mass Effects"
C       By: S. Kretzer, H.L. Lai, F. Olness, W.K. Tung
C       Phys. Rev. D69:114005(2004), hep-ph/0307022
C
C   Ref[5]: "Heavy Quark Mass Effects in Deep Inelastic Scattering and Global QCD Analysis"
C       By : W.K. Tung, H.L. Lai, A. Belyaev, J. Pumplin, D. Stump, C.-P. Yuan
C       JHEP 0702:053(2007), hep-ph/0611254
C
C   Ref[6]: "The Strange Parton Distribution of Nucleon: Global Analysis and Applications"
C       By : H.L. Lai, P. Nadolsky, J. Pumplin, D. Stump, W.K. Tung, C.-P. Yuan
C       JHEP 0704:089,2007, hep-ph/0702268
C
C   Ref[7]: "The Charm Content of the Nucleon"
C       By : J. Pumplin, H.L. Lai, W.K. Tung
C       Phys.Rev.D75:054029,2007, hep-ph/0701220

C   Ref[8]: "Implications of CTEQ global analysis for collider observables"
C       By : P. M. Nadolsky, H.-L. Lai, Q.-H. Cao, J. Huston, J. Pumplin, D. R. Stump, W.-K. Tung, C.-P. Yuan
C       arXiv:0802.0007 [hep-ph], submitted to Phys. Rev. D.
C

C   This package contains
C   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
C   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
C   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].
C   (4) 5 special sets for strangeness study from Ref[3].
C   (5) 1 special set for heavy quark study from Ref[4].
C   (6) CTEQ6.5M and its 40 up/down eigenvector sets from Ref[5].
C   (7) 8 sets of PDFs resulting from the strangeness study, Ref[6].
C   (8) 7 sets of PDFs resulting from the charm study, Ref[7].
C   (9) CTEQ6.6M and its 44 up/down eigenvector sets from Ref[8].
C  (10) Fits with nonperturbative charm from the study in  Ref[8].
C  (11) Fits with alternative values of the strong coupling strength from the study in Ref[8].


C  Details about the calling convention are:
C --------------------------------------------------------------------------------
C  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File   Ref
C ================================================================================
C Standard, "best-fit", sets:
C --------------------------
C   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl    [1]
C   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl    [1]
C   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl    [1]
C   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl   [1]
C 200    CTEQ6.1M: updated CTEQ6M (see below, under "uncertainty" section)     [2]
C 400    CTEQ6.6M; the 2008 set (see below, under "uncertainty" section)       [8]
C
C --------------------------
C  Special sets with nonperturbative charm at Q_0=1.3 GeV from Ref [8]
C --------------------------
C 450    CTEQ6.6C1   BHPS model for IC    0.118     326   226    ctq66.c1.pds
C 451    CTEQ6.6C2   BHPS model for IC    0.118     326   226    ctq66.c2.pds
C 452    CTEQ6.6C3   Sea-like model       0.118     326   226    ctq66.c3.pds
C 453    CTEQ6.6C4   Sea-like model       0.118     326   226    ctq66.c4.pds
C     Momentum Fraction carried by c+cbar=2c at Q0=1.3 GeV:
C    Iset:     451  452   453   454
C Mom. frac:  0.01 0.035  0.01  0.035


C --------------------------
C  Special CTEQ6.6 sets with alternative values of strong coupling strength [8]
C --------------------------
C 460    CTEQ6.6A1                        0.125           328    ctq66.a1.pds
C 461    CTEQ6.6A2                        0.122           281    ctq66.a2.pds
C 462    CTEQ6.6A3                        0.114           179    ctq66.a3.pds
C 463    CTEQ6.6A4                        0.112           159    ctq66.a4.pds

C --------------------------
C Special sets for strangeness study:  Ref.[3]
C --------------------------
C  11    CTEQ6A   Class A                 0.118     326   226    cteq6sa.pds
C  12    CTEQ6B   Class B                 0.118     326   226    cteq6sb.pds
C  13    CTEQ6C   Class C                 0.118     326   226    cteq6sc.pds
C  14    CTEQ6B+  Large [S-]              0.118     326   226    cteq6sb+.pds
C  15    CTEQ6B-  Negative [S-]           0.118     326   226    cteq6sb-.pds
C --------------------------
C Special set for Heavy Quark study:   Ref.[4]
C --------------------------
C  21    CTEQ6HQ                          0.118     326   226    cteq6hq.pds
C --------------------------
C Released sets for strangeness study:  Ref.[6]
C -------------------------- s=sbr
C  30    CTEQ6.5S0   Best-fit             0.118     326   226    ctq65.s+0.pds
C  31    CTEQ6.5S1   Low s+               0.118     326   226    ctq65.s+1.pds
C  32    CTEQ6.5S2   High s+              0.118     326   226    ctq65.s+2.pds
C  33    CTEQ6.5S3   Alt Low s+           0.118     326   226    ctq65.s+3.pds
C  34    CTEQ6.5S4   Alt High s+          0.118     326   226    ctq65.s+4.pds
C -------------------------- s!=sbr
C          strangeness asymmetry <x>_s-
C  35    CTEQ6.5S-0  Best-fit    0.0014    0.118     326   226    ctq65.s-0.pds
C  36    CTEQ6.5S-1  Low        -0.0010    0.118     326   226    ctq65.s-1.pds
C  37    CTEQ6.5S-2  High        0.0050    0.118     326   226    ctq65.s-2.pds
C --------------------------
C Released sets for charm study:  Ref.[7]
C --------------------------
C  40    CTEQ6.5C0   no intrinsic charm   0.118     326   226    ctq65.c0.pds
C  41    CTEQ6.5C1   BHPS model for IC    0.118     326   226    ctq65.c1.pds
C  42    CTEQ6.5C2   BHPS model for IC    0.118     326   226    ctq65.c2.pds
C  43    CTEQ6.5C3   Meson cloud model    0.118     326   226    ctq65.c3.pds
C  44    CTEQ6.5C4   Meson cloud model    0.118     326   226    ctq65.c4.pds
C  45    CTEQ6.5C5   Sea-like model       0.118     326   226    ctq65.c5.pds
C  46    CTEQ6.5C6   Sea-like model       0.118     326   226    ctq65.c6.pds
C
C     Momentum Fraction carried by c,cbar at Q0=1.3 GeV:
C    Iset:charm  ,cbar     | Iset:charm  ,cbar     | Iset:charm  ,cbar
C    41: 0.002857,0.002857 | 43: 0.003755,0.004817 | 45: 0.005714,0.005714
C    42: 0.010000,0.010000 | 44: 0.007259,0.009312 | 46: 0.012285,0.012285
C
C ============================================================================
C For uncertainty calculations using eigenvectors of the Hessian:
C ---------------------------------------------------------------
C     central + 40 up/down sets along 20 eigenvector directions
C                             -----------------------------
C                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
C                             -----------------------
C  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
C             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
C                              -----------------------
C  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
C             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects, Ref[5]:  central fit: CTEQ6.5M (=CTEQ65.00)
C                              -----------------------
C  3xx  CTEQ65.xx  +/- sets               0.118     326   226    ctq65.xx.pds
C        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 300      is CTEQ65.00 (=CTEQ6.5M),
C             301/302 are CTEQ65.01/02, +/- sets of 1st eigenvector, ... etc.
C        ====================================================================
C                Version with mass effects and free strangeness, Ref[8]:
C                central fit: CTEQ6.6M (=CTEQ66.00)
C                              -----------------------
C  4xx  CTEQ66.xx  +/- sets               0.118     326   226    ctq66.xx.pds
C        where xx = 01-44: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
C        e.g. 400      is CTEQ66.00 (=CTEQ6.6M),
C             401/402 are CTEQ66.01/02, +/- sets of 1st eigenvector, ... etc.

C ===========================================================================
C   ** ALL fits are obtained by using the same coupling strength
C   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
C   which uses the LO running \alpha_s and its value determined from the fit.
C   For the LO fits, the evolution of the PDF and the hard cross sections are
C   calculated at LO.  More detailed discussions are given in the references.
C
C   The table grids are generated for
C    *  10^-8 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.6 series;
C    *  10^-7 < x < 1 and 1.3 < Q < 10^5 (GeV) for CTEQ6.5S/C series;
C    *  10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV) for CTEQ6, CTEQ6.1 series;
C
C   PDF values outside of the above range are returned using extrapolation.
C   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
C   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
C   which is defined as the bottom quark mass, whenever it can be applied.
C
C   The Table_Files are assumed to be in the working directory.
C
C   Before using the PDF, it is necessary to do the initialization by
C       Call SetCtq6(Iset)
C   where Iset is the desired PDF specified in the above table.
C
C   The function Ctq6Pdf (Iparton, X, Q)
C   returns the parton distribution inside the proton for parton [Iparton]
C   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
C   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
C                            for (b, c, s, d, u, g, u_bar, ..., b_bar),
C
C   For detailed information on the parameters used, e.q. quark masses,
C   QCD Lambda, ... etc.,  see info lines at the beginning of the
C   Table_Files.
C
C   These programs, as provided, are in double precision.  By removing the
C   "Implicit Double Precision" lines, they can also be run in single
C   precision.
C
C   If you have detailed questions concerning these CTEQ6 distributions,
C   or if you find problems/bugs using this package, direct inquires to
C   nadolsky@pa.msu.edu, pumplin@pa.msu.edu or tung@pa.msu.edu.
C
C===========================================================================

      Function Ctq6Pdf(Iparton, X, Q)
csp      Implicit Double Precision (A-H,O-Z)
      Logical Warn
      Common
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / QCDtable /  Alambda, Nfl, Iorder

      Data Warn /.true./
      save Warn

      If (X .lt. 0d0 .or. X .gt. 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
        Ctq6Pdf = 0D0
        Return
      Endif

      If (Q .lt. Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
        Stop
      Endif

      If ((Iparton .lt. -NfMx .or. Iparton .gt. NfMx)) Then
         If (Warn) Then
C        put a warning for calling extra flavor.
             Warn = .false.
             Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
             Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
         Endif
         Ctq6Pdf = 0D0
         Return
      Endif

      Ctq6Pdf = PartonX6 (Iparton, X, Q)
      if (Ctq6Pdf.lt.0.D0) Ctq6Pdf = 0.D0

      Return

C                             ********************
      End

      Subroutine SetCtq6 (Iset)
      Implicit Double Precision (A-H,O-Z)
      integer mxho
      parameter(mxho=10)
      character*500 fnch,fnhi,fndt,fnhm,fnii,fnid,fnie,fnrj,fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnio,fnho,fn3g,fn3p,fn3d
     *,fn3f,fn3f1,fn3f2,fn3f3,fn3f4,fn3f5,fnhpf
      common/fname/ fnch,fnhi,fndt,fnhm,fnii,fnid,fnie,fnrj,fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnio,fnho(mxho),fn3g,fn3p,fn3d
     *,fn3f,fn3f1,fn3f2,fn3f3,fn3f4,fn3f5,fnhpf
      integer nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnio,nfnho,nfn3g,nfnhm
     *,nfn3p,nfn3d,nfn3f,nfn3f1,nfn3f2,nfn3f3,nfn3f4,nfn3f5,nfnhpf
      common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnio,nfnho(mxho),nfn3g,nfnhm
     *,nfn3p,nfn3d,nfn3f,nfn3f1,nfn3f2,nfn3f3,nfn3f4,nfn3f5,nfnhpf
      Parameter (Isetmax0=8)
      Character Flnm(Isetmax0)*6, nn*3, Tablefile*540
      Logical fmtpds
      Data (Flnm(I), I=1,Isetmax0)
     > / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.','cteq6s'
     >  ,'ctq65.', 'ctq66.' /
      Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
      Data Isetmin2,Isetmax2 /200,240/
      Data Isetmin3,Isetmax3 /300,340/
      Data Isetmin4,Isetmax4 /400,444/
      Data IsetminS,IsetmaxS /11,15/
      Data IsetmnSp07,IsetmxSp07 /30,34/
      Data IsetmnSm07,IsetmxSm07 /35,37/
      Data IsetmnC07,IsetmxC07 /40,46/
      Data IsetmnC08,IsetmxC08 /450,453/
      Data IsetmnAS08,IsetmxAS08 /460,463/

      Data IsetHQ /21/
      Common /Setchange6/ Isetch
      save

C             If data file not initialized, do so.
      If(Iset.ne.Isetold) then
        fmtpds=.true.

        If (Iset.ge.Isetmin0 .and. Iset.le.3) Then
C                                                  Iset = 1,2,3 for 6m, 6d, 6l
          fmtpds=.false.
          Tablefile=fngrv(1:nfngrv)//Flnm(Iset)//'.tbl'
        Elseif (Iset.eq.4) Then
C                                                             4  (2nd LO fit)
          fmtpds=.false.
          Tablefile=fngrv(1:nfngrv)//Flnm(Iset)//'1.tbl'
        Elseif (Iset.ge.Isetmin1 .and. Iset.le.Isetmax1) Then
C                                                               101 - 140
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(1)//nn//'.tbl'
        Elseif (Iset.ge.Isetmin2 .and. Iset.le.Isetmax2) Then
C                                                               200 - 240
          fmtpds=.false.
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(5)//nn(2:3)//'.tbl'
        Elseif (Iset.ge.IsetminS .and. Iset.le.IsetmaxS) Then
C                                                               11 - 15
          If(Iset.eq.11) then
            Tablefile=fngrv(1:nfngrv)//Flnm(6)//'a.pds'
          Elseif(Iset.eq.12) then
            Tablefile=fngrv(1:nfngrv)//Flnm(6)//'b.pds'
          Elseif(Iset.eq.13) then
            Tablefile=fngrv(1:nfngrv)//Flnm(6)//'c.pds'
          Elseif(Iset.eq.14) then
            Tablefile=fngrv(1:nfngrv)//Flnm(6)//'b+.pds'
          Elseif(Iset.eq.15) then
            Tablefile=fngrv(1:nfngrv)//Flnm(6)//'b-.pds'
          Endif
        Elseif (Iset.eq.IsetHQ) Then
C                                                               21
          TableFile=fngrv(1:nfngrv)//'cteq6hq.pds'
        Elseif (Iset.ge.IsetmnSp07 .and. Iset.le.IsetmxSp07) Then
C                                                    (Cteq6.5S)  30 - 34
          write(nn,'(I2)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(7)//'s+'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnSm07 .and. Iset.le.IsetmxSm07) Then
C                                                    (Cteq6.5S)  35 - 37
          Is = Iset - 5
          write(nn,'(I2)') Is
          Tablefile=fngrv(1:nfngrv)//Flnm(7)//'s-'//nn(2:2)//'.pds'
        Elseif (Iset.ge.IsetmnC07 .and. Iset.le.IsetmxC07) Then
C                                                    (Cteq6.5C)  40 - 46
          write(nn,'(I2)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(7)//'c'//nn(2:2)//'.pds'
        Elseif (Iset.ge.Isetmin3 .and. Iset.le.Isetmax3) Then
C                                                    (Cteq6.5)  300 - 340
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(7)//nn(2:3)//'.pds'
        Elseif (Iset.ge.Isetmin4 .and. Iset.le.Isetmax4) Then
C                                                    (Cteq6.6)  400 - 444
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(8)//nn(2:3)//'.pds'
        Elseif (Iset.ge.IsetmnC08 .and. Iset.le.IsetmxC08) Then
C                                                   (Cteq6.6C)  450 - 453
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(8)//'c'//nn(3:3)//'.pds'
        Elseif (Iset.ge.IsetmnAS08 .and. Iset.le.IsetmxAS08) Then
C                                                   (Cteq6.6AS)  460 - 463
          write(nn,'(I3)') Iset
          Tablefile=fngrv(1:nfngrv)//Flnm(8)//'a'//nn(3:3)//'.pds'
        Else
          Print *, 'Invalid Iset number in SetCtq6 :', Iset
          Stop
        Endif
        IU= NextUn()
        Open(IU, File=Tablefile, Status='OLD', Err=100)
        Call Readpds (IU,fmtpds)
        Close (IU)
        Isetold=Iset
        Isetch=1
      Endif
      Return

 100  Print *, ' Data file ', Tablefile
     >       , ' cannot be opened in SetCtq6!!'
      Stop
C                             ********************
      End

      Subroutine Readpds (Nu,fmtpds)
cs[      Implicit Double Precision (A-H,O-Z)
      Character Line*80
      Logical fmtpds
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > / QCDtable /  Alambda, Nfl, Iorder
     > / Masstbl / Amass(6)

      Read  (Nu, '(A)') Line
      Read  (Nu, '(A)') Line
      Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
      Iorder = Nint(Dr)
      Nfl = Nint(Fl)
      Alambda = Al

      Read  (Nu, '(A)') Line
      If(fmtpds) then
C                                               This is the .pds (WKT) format
        Read  (Nu, *) N0, N0, N0, NfMx, MxVal, N0
        If(MxVal.gt.MaxVal) MxVal=3 !old .pds format (read in KF, not MxVal)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) NX,  NT, N0, NG, N0

        Read  (Nu, '(A)') (Line,I=1,NG+2)
        Read  (Nu, *) QINI, QMAX, (aa,TV(I), I =0, NT)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
        XV(0)=0D0
      Else
C                                               This is the old .tbl (HLL) format
         MxVal=2
         Read  (Nu, *) NX,  NT, NfMx

         Read  (Nu, '(A)') Line
         Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

         Read  (Nu, '(A)') Line
         Read  (Nu, *) XMIN, (XV(I), I =0, NX)

         Do 11 Iq = 0, NT
            TV(Iq) = Log(Log (TV(Iq) /Al))
 11      Continue
      Endif

      Nblk = (NX+1) * (NT+1)
      Npts =  Nblk  * (NfMx+1+MxVal)
      Read  (Nu, '(A)') Line
      Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

      Return
C                        ****************************
      End

      Function PartonX6 (IPRTN, XX, QQ)

c  Given the parton distribution function in the array U in
c  COMMON / PEVLDT / , this routine interpolates to find
c  the parton distribution at an arbitray point in x and q.
c
csp      Implicit Double Precision (A-H,O-Z)
csp
      PARAMETER (MXX = 201, MXQ = 25, MXF = 6, MaxVal=4)
      PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

      Common
     > / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX)
     > / CtqPar2 / Nx, Nt, NfMx, MxVal
     > / XQrange / Qini, Qmax, Xmin
     > /Setchange6/ Isetch

      Dimension fvec(4), fij(4)
      Dimension xvpow(0:mxx)
      Data OneP / 1.00001 /
      Data xpow / 0.3d0 /       !**** choice of interpolation variable
      Data nqvec / 4 /
      Data ientry / 0 /
      Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
      Save xvpow
      Save X, Q, JX, JQ, JLX, JLQ
      Save ss, const1, const2, const3, const4, const5, const6
      Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
      Save tmp1, tmp2, tdet

      If((XX.eq.X).and.(QQ.eq.Q)) goto 99
c store the powers used for interpolation on first call...
      if(Isetch .eq. 1) then
         Isetch = 0

         xvpow(0) = 0D0
         do i = 1, nx
            xvpow(i) = xv(i)**xpow
         enddo
      endif

      X = XX
      Q = QQ
      tt = log(log(Q/Al))

c      -------------    find lower end of interval containing x, i.e.,
c                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
      JLx = -1
      JU = Nx+1
 11   If (JU-JLx .GT. 1) Then
         JM = (JU+JLx) / 2
         If (X .Ge. XV(JM)) Then
            JLx = JM
         Else
            JU = JM
         Endif
         Goto 11
      Endif
C                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
C                           |---|---|---|...|---|-x-|---|...|---|---|
C                     x     0  Xmin               x                 1
C
      If     (JLx .LE. -1) Then
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
        Stop
      ElseIf (JLx .Eq. 0) Then
         Jx = 0
      Elseif (JLx .LE. Nx-2) Then

C                For interrior points, keep x in the middle, as shown above
         Jx = JLx - 1
      Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then

C                  We tolerate a slight over-shoot of one (OneP=1.00001),
C              perhaps due to roundoff or whatever, but not more than that.
C                                      Keep at least 4 points >= Jx
         Jx = JLx - 2
      Else
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
        Stop
      Endif
C          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

C                       This is the variable to be interpolated in
      ss = x**xpow

      If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then

c     initiation work for "interior bins": store the lattice points in s...
      svec1 = xvpow(jx)
      svec2 = xvpow(jx+1)
      svec3 = xvpow(jx+2)
      svec4 = xvpow(jx+3)

      s12 = svec1 - svec2
      s13 = svec1 - svec3
      s23 = svec2 - svec3
      s24 = svec2 - svec4
      s34 = svec3 - svec4

      sy2 = ss - svec2
      sy3 = ss - svec3

c constants needed for interpolating in s at fixed t lattice points...
      const1 = s13/s23
      const2 = s12/s23
      const3 = s34/s23
      const4 = s24/s23
      s1213 = s12 + s13
      s2434 = s24 + s34
      sdet = s12*s34 - s1213*s2434
      tmp = sy2*sy3/sdet
      const5 = (s34*sy2-s2434*sy3)*tmp/s12
      const6 = (s1213*sy2-s12*sy3)*tmp/s34

      EndIf

c         --------------Now find lower end of interval containing Q, i.e.,
c                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
      JLq = -1
      JU = NT+1
 12   If (JU-JLq .GT. 1) Then
         JM = (JU+JLq) / 2
         If (tt .GE. TV(JM)) Then
            JLq = JM
         Else
            JU = JM
         Endif
         Goto 12
       Endif

      If     (JLq .LE. 0) Then
         Jq = 0
      Elseif (JLq .LE. Nt-2) Then
C                                  keep q in the middle, as shown above
         Jq = JLq - 1
      Else
C                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

      Endif
C                                   This is the interpolation variable in Q

      If (JLq.GE.1 .and. JLq.LE.Nt-2) Then
c                                        store the lattice points in t...
      tvec1 = Tv(jq)
      tvec2 = Tv(jq+1)
      tvec3 = Tv(jq+2)
      tvec4 = Tv(jq+3)

      t12 = tvec1 - tvec2
      t13 = tvec1 - tvec3
      t23 = tvec2 - tvec3
      t24 = tvec2 - tvec4
      t34 = tvec3 - tvec4

      ty2 = tt - tvec2
      ty3 = tt - tvec3

      tmp1 = t12 + t13
      tmp2 = t24 + t34

      tdet = t12*t34 - tmp1*tmp2

      EndIf


c get the pdf function values at the lattice points...

 99   If (Iprtn .Gt. MxVal) Then
         Ip = - Iprtn
      Else
         Ip = Iprtn
      EndIf
      jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

      Do it = 1, nqvec

         J1  = jtmp + it*(NX+1)

       If (Jx .Eq. 0) Then
C                          For the first 4 x points, interpolate x^2*f(x,Q)
C                           This applies to the two lowest bins JLx = 0, 1
C            We can not put the JLx.eq.1 bin into the "interrior" section
C                           (as we do for q), since Upd(J1) is undefined.
         fij(1) = 0
         fij(2) = Upd(J1+1) * XV(1)**2
         fij(3) = Upd(J1+2) * XV(2)**2
         fij(4) = Upd(J1+3) * XV(3)**2
C
C                 Use Polint which allows x to be anywhere w.r.t. the grid

         Call Polint4F (XVpow(0), Fij(1), ss, Fx)

         If (x .GT. 1D-15)  Fvec(it) =  Fx / x**2
C                                              Pdf is undefined for x.eq.0
       ElseIf  (JLx .Eq. Nx-1) Then
C                                                This is the highest x bin:

        Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

        Fvec(it) = Fx

       Else
C                       for all interior points, use Jon's in-line function
C                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
         sf2 = Upd(J1+1)
         sf3 = Upd(J1+2)

         g1 =  sf2*const1 - sf3*const2
         g4 = -sf2*const3 + sf3*const4

         Fvec(it) = (const5*(Upd(J1)-g1)
     &               + const6*(Upd(J1+3)-g4)
     &               + sf2*sy3 - sf3*sy2) / s23

       Endif

      enddo
C                                   We now have the four values Fvec(1:4)
c     interpolate in t...

      If (JLq .LE. 0) Then
C                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

      ElseIf (JLq .GE. Nt-1) Then
C                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
      Else
C                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
C       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
C                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12
     &    +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
      EndIf

      PartonX6 = ff

      Return
C                                       ********************
      End

      SUBROUTINE POLINT4F (XA,YA,X,Y)

csp      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
C  but assuming N=4, and ignoring the error estimation.
C  suggested by Z. Sullivan.
      DIMENSION XA(*),YA(*)

      H1=XA(1)-X
      H2=XA(2)-X
      H3=XA(3)-X
      H4=XA(4)-X

      W=YA(2)-YA(1)
      DEN=W/(H1-H2)
      D1=H2*DEN
      C1=H1*DEN

      W=YA(3)-YA(2)
      DEN=W/(H2-H3)
      D2=H3*DEN
      C2=H2*DEN

      W=YA(4)-YA(3)
      DEN=W/(H3-H4)
      D3=H4*DEN
      C3=H3*DEN

      W=C2-D1
      DEN=W/(H1-H3)
      CD1=H3*DEN
      CC1=H1*DEN

      W=C3-D2
      DEN=W/(H2-H4)
      CD2=H4*DEN
      CC2=H2*DEN

      W=CC2-CD1
      DEN=W/(H1-H4)
      DD1=H4*DEN
      DC1=H1*DEN

      If((H3+H4).lt.0D0) Then
         Y=YA(4)+D3+CD2+DD1
      Elseif((H2+H3).lt.0D0) Then
         Y=YA(3)+D2+CD1+DC1
      Elseif((H1+H2).lt.0D0) Then
         Y=YA(2)+C2+CD1+DC1
      ELSE
         Y=YA(1)+C1+CC1+DC1
      ENDIF

      RETURN
C               *************************
      END

      Function NextUn()
C                                 Returns an unallocated FORTRAN i/o unit.
      Logical EX
C
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUn = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '
C               *************************
      End
C

csp------------------------------------------------------------------
      function cteq6(x,qqs,qq,icq,iq)
csp------------------------------------------------------------------
csp PDF from CTEQ6
csp------------------------------------------------------------------
       common/cpicteq6/npicteq6
       data npicteq6/0/

       dum=qq
       idum=icq

       npicteq6=npicteq6+1
c       if(npicteq6.eq.1)call SetCtq6(400)     !Cteq 6.6
       if(npicteq6.eq.1)call SetCtq6(3)    !Cteq 6.1

       if(x.gt..999999.or.qqs.le.0.25)then
        cteq6=0.
        return
       endif

       if(iq.eq.0)then
       cteq6=x*Ctq6Pdf(0,x,qqs)
       elseif(iq.eq.1)then
       cteq6=x*max(0.,Ctq6Pdf(1,x,qqs)-Ctq6Pdf(-1,x,qqs))
       elseif(iq.eq.2)then
       cteq6=x*max(0.,Ctq6Pdf(2,x,qqs)-Ctq6Pdf(-2,x,qqs))
       elseif(iq.eq.3)then
       cteq6=x*max(0.,Ctq6Pdf(3,x,qqs)-Ctq6Pdf(-3,x,qqs))
       elseif(iq.eq.4)then                                  !bg
       cteq6=x*max(0.,Ctq6Pdf(4,x,qqs))   !tp
       elseif(iq.eq.5)then                                  !bg
       cteq6=x*max(0.,Ctq6Pdf(5,x,qqs))   !tp
       elseif(iq.eq.-4)then                                 !bg
       cteq6=x*Ctq6Pdf(-4,x,qqs)                            !bg
       elseif(iq.eq.-5)then                                 !tp
       cteq6=x*Ctq6Pdf(-5,x,qqs)                            !tp
       elseif(iq.eq.-3)then
       cteq6=x*Ctq6Pdf(-3,x,qqs)
       elseif(iq.eq.-1)then
       cteq6=x*Ctq6Pdf(-1,x,qqs)
       elseif(iq.eq.-2)then
       cteq6=x*Ctq6Pdf(-2,x,qqs)
       else
       cteq6=0.
       endif

      end
c-----------------------------------------------------------------------
      double precision function psjcteq6x(t,qt,xx1,xx2,s)
c-----------------------------------------------------------------------
ctp 20171005 Exactly the same results using psbori directly or same way
c of counting quarks than in ffsig : Formula OK
c-----------------------------------------------------------------------
#include "sem.h"
      double precision ffborn,t,s!,g1,ub1,u1,db1,d1,sb1,s1
c     *                           ,g2,ub2,u2,db2,d2,sb2,s2,psbori

      g1=cteq6(xx1,qt,0.,2,0)
      ub1=cteq6(xx1,qt,0.,2,-1)      !aq sea
      uv1=cteq6(xx1,qt,0.,2,1)       !q val
      db1=cteq6(xx1,qt,0.,2,-2)      !aq sea
      dv1=cteq6(xx1,qt,0.,2,2)        !q val 
      sb1=cteq6(xx1,qt,0.,2,-3)      !aq sea
      sea1=(sb1+ub1+db1)/3.     !q sea
      ch1=cteq6(xx1,qt,0.,2,-4)
      bt1=cteq6(xx1,qt,0.,2,-5)
      g2=cteq6(xx2,qt,0.,2,0)
      ub2=cteq6(xx2,qt,0.,2,-1)
      uv2=cteq6(xx2,qt,0.,2,1)
      db2=cteq6(xx2,qt,0.,2,-2)
      dv2=cteq6(xx2,qt,0.,2,2)
      sb2=cteq6(xx2,qt,0.,2,-3)
      sea2=(sb2+ub2+db2)/3.
      ch2=cteq6(xx2,qt,0.,2,-4)
      bt2=cteq6(xx2,qt,0.,2,-5)


      call ffvalues(qt, g1,uv1,dv1,sea1,g2,uv2,dv2,sea2 ! <--in
     .      ,ch1,bt1,ch2,bt2
     .      ,gg,gq,qq,qa,qqp,gc,gb,qc,qb,cc,cac,bb,bab,cb)! <--out
      psjcteq6x=ffborn(0,99,99d0,s,t,gg,gq,qq,qa,qqp ,1,3, 3 ) 
     *         +ffborn(0,99,99d0,s,t,0.,gc,0.,0. ,qc ,1,3, 4 )
     *         +ffborn(0,99,99d0,s,t,0.,gb,0.,0. ,qb ,1,3, 5 )
     *         +ffborn(0,99,99d0,s,t,0.,0.,cc,cac,0. ,1,3, 6 )
     *         +ffborn(0,99,99d0,s,t,0.,0.,bb,bab,cb ,1,3, 7 )

c      print *,'psj',s,t,qt,xx1,xx2,g1,uv1,dv1,sea1,g2,uv2,dv2,sea2,gc,qc
c     *       ,psjcteq6x


      return
      end

cc-----------------------------------------------------------------------
c      double precision function psjcteq6y(t,qt,xx1,xx2,s)
cc-----------------------------------------------------------------------
cctp 20171005 Exactly the same results using psbori directly or same way
cc of counting quarks than in ffsig : Formula OK
cc-----------------------------------------------------------------------
c#include "sem.h"
c      double precision ffborn,t,s!,g1,ub1,u1,db1,d1,sb1,s1
cc     *                           ,g2,ub2,u2,db2,d2,sb2,s2,psbori
c
c      g1=cteq6(xx1,qt,0.,2,0)
c      ub1=cteq6(xx1,qt,0.,2,-1)      !aq sea
c      u1=cteq6(xx1,qt,0.,2,1)+ub1    !q val + q sea
c      db1=cteq6(xx1,qt,0.,2,-2)      !aq sea
c      d1=cteq6(xx1,qt,0.,2,2)+db1   !q val + q sea
c      sb1=cteq6(xx1,qt,0.,2,-3)      !aq sea
c      s1=sb1                         !q sea
c      g2=cteq6(xx2,qt,0.,2,0)
c      ub2=cteq6(xx2,qt,0.,2,-1)
c      u2=cteq6(xx2,qt,0.,2,1)+ub2
c      db2=cteq6(xx2,qt,0.,2,-2)
c      d2=cteq6(xx2,qt,0.,2,2)+db2
c      sb2=cteq6(xx2,qt,0.,2,-3)
c      s2=sb2
c
cc      g1=dble(cteq6(xx1,qt,0.,2,0))
cc      ub1=dble(cteq6(xx1,qt,0.,2,-1))
cc      u1=dble(cteq6(xx1,qt,0.,2,1))+ub1
cc      db1=dble(cteq6(xx1,qt,0.,2,-2))
cc      d1=dble(cteq6(xx1,qt,0.,2,2))+db1
cc      sb1=dble(cteq6(xx1,qt,0.,2,-3))
cc      s1=sb1
cc      g2=dble(cteq6(xx2,qt,0.,2,0))
cc      ub2=dble(cteq6(xx2,qt,0.,2,-1))
cc      u2=dble(cteq6(xx2,qt,0.,2,1))+ub2
cc      db2=dble(cteq6(xx2,qt,0.,2,-2))
cc      d2=dble(cteq6(xx2,qt,0.,2,2))+db2
cc      sb2=dble(cteq6(xx2,qt,0.,2,-3))
cc      s2=sb2
c
c
c       psjcteq6x=ffborn(0,99,99d0,s,t,  g1*g2                   !gg
c
c     *,(g2*(u1+ub1+d1+db1+s1+sb1)+g1*(u2+ub2+d2+db2+s2+sb2))   !gq
c
c     *,(u1*u2+ub1*ub2+d1*d2+db1*db2+s1*s2+sb1*sb2)            !qq
c
c     *,(u1*ub2+ub1*u2+d1*db2+db1*d2+s1*sb2+sb1*s2)            !qa
c
c     *,((u1+ub1)*(d2+db2+s2+sb2)+(u2+ub2)*(d1+db1+s1+sb1)+
c     *(d1+db1)*(u2+ub2+s2+sb2)+(d2+db2)*(u1+ub1+s1+sb1)+
c     *(s1+sb1)*(u2+ub2+d2+db2)+(s2+sb2)*(u1+ub1+d1+db1))      !qqp
c     *,1,3,3)
c
cc      uv1=cteq6(xx1,qt,0.,2,1)       !q val
cc      dv1=cteq6(xx1,qt,0.,2,2)        !q val 
cc      sea1=(sb1+ub1+db1)/3.           !q sea
cc      uv2=cteq6(xx2,qt,0.,2,1)
cc      dv2=cteq6(xx2,qt,0.,2,2)
cc      sea2=(sb2+ub2+db2)/3.
cc       print *,'pdf',s,t,xx1,xx2,qt                 !gg
cc     *,(g2*(u1+ub1+d1+db1+s1+sb1)+g1*(u2+ub2+d2+db2+s2+sb2))   !gq
cc     *,(uv1+dv1+6.*sea1)*g2 + g1*(uv2+dv2+6.*sea2) 
cc     *,uv1,dv1,sea1,g2,g1,uv2,dv2,sea2
c
c       
cc      psjcteq6x=g1*g2*(psbori(99,s,t,0,0,1)+psbori(99,s,s-t,0,0,1)
cc     *+psbori(99,s,t,0,0,2)+psbori(99,s,s-t,0,0,2))/2.
cc
cc     *+(psbori(99,s,t,0,1,1)+psbori(99,s,s-t,0,1,1))*
cc     *(g2*(u1+ub1+d1+db1+s1+sb1)+g1*(u2+ub2+d2+db2+s2+sb2))
cc
cc     *+(psbori(99,s,t,1,1,1)+psbori(99,s,s-t,1,1,1))/2.*
cc     *(u1*u2+ub1*ub2+d1*d2+db1*db2+s1*s2+sb1*sb2)
cc
cc     *+(psbori(99,s,t,1,-1,1)+psbori(99,s,s-t,1,-1,1)
cc     *+psbori(99,s,t,1,-1,2)+
cc     *psbori(99,s,s-t,1,-1,2)+psbori(99,s,t,1,-1,3)
cc     *+psbori(99,s,s-t,1,-1,3))*
cc     *(u1*ub2+ub1*u2+d1*db2+db1*d2+s1*sb2+sb1*s2)
cc
cc     *+(psbori(99,s,t,1,2,1)+psbori(99,s,s-t,1,2,1))*
cc     *((u1+ub1)*(d2+db2+s2+sb2)+(u2+ub2)*(d1+db1+s1+sb1)+
cc     *(d1+db1)*(u2+ub2+s2+sb2)+(d2+db2)*(u1+ub1+s1+sb1)+
cc     *(s1+sb1)*(u2+ub2+d2+db2)+(s2+sb2)*(u1+ub1+d1+db1))
c
c      return
c      end

c-----------------------------------------------------------------------bg
      subroutine Cparton        !bg
c--------------------------------------------------------------------------
c     calculation for charm production based on parton model
c     divide by 2 because psbori(99,s,t,...)+psbori(99,s,s-t=u,...)?
c-------------------------------------------------------------------------
#include "sem.h"
#include "aaa.h"
#include "ems.h"
      double precision ggcteq6x,cgcteq6x,gglcteq6x,qaqteq6x,qpcteq6x
      external ggcteq6x,cgcteq6x,gglcteq6x,qaqteq6x,qpcteq6x
      dimension sgg(50),scg(50),sggl(50),sqaq(50),sqpc(50),ssum(50)

      sx=engy**2
      y2=1.
      do i=1,40
        pt=real(i)
        sgg(i)=psjpdf(ggcteq6x,pt**2,sx,y2)/2.
        sggl(i)=psjpdf(gglcteq6x,pt**2,sx,y2)/2.
        scg(i)=psjpdf(cgcteq6x,pt**2,sx,y2)/2.
        sqaq(i)=psjpdf(qaqteq6x,pt**2,sx,y2)/2.
        sqpc(i)=psjpdf(qpcteq6x,pt**2,sx,y2)/2.
        ssum(i)=sgg(i)+scg(i)+sqaq(i)+sqpc(i)/2.
      enddo

      open(76,file='charm-parton.histo')
      write(76,'(a)')       'orientation landscape'
      write(76,'(a)')       'newpage'
      write(76,'(a)')       'zone 1 1 1'
      write(76,'(a)')'set scale 2 set scalel 0.5  set scalem 0.9'
      write(76,'(a)')       'set basl 0.01'
      write(76,'(a)')        'cd .'
      write(76,'(a)')  'openhisto name gg'
      write(76,'(a)')  'htyp lgu xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),sgg(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0-'

      write(76,'(a)')  'openhisto name ggl'
      write(76,'(a)')  'htyp lbu xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),sggl(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0-'

      write(76,'(a)')  'openhisto name cg'
      write(76,'(a)')  'htyp lru xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),scg(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0-'

      write(76,'(a)')  'openhisto name qpc'
      write(76,'(a)')  'htyp lxu xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),sqpc(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0-'

      write(76,'(a)')  'openhisto name qaq'
      write(76,'(a)')  'htyp lyu xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),sqaq(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0-'

      write(76,'(a)')  'openhisto name sum'
      write(76,'(a)')  'htyp lku xmod lin ymod log'
      write(76,'(a)')  'xrange 0.0 70.0 yrange auto auto'
      write(76,'(a)')  'text 0 0 ""title""'
      write(76,'(a)')  'txt "xaxis pt "'
      write(76,'(a)')  'txt "yaxis dn/dpt  "'
      write(76,'(a)')  'histoweight 1'
      write(76,'(a)')  ' array 2'
      do i=1,40
        write(76,*) real(i),ssum(i)*factk/sigine*10*.0389/y2
      enddo
      write(76,'(a)')  'endarray'
      write(76,'(a)')  'closehisto plot 0'

      close(76)

      end

c-----------------------------------------------------------------------bg
      double precision function ggcteq6x(t,qt,xx1,xx2,s) !bg
c-----------------------------------------------------------------------------
c     gg->cc~
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s

      g1=cteq6(xx1,qt,0.,2,0)
      g2=cteq6(xx2,qt,0.,2,0)

      ggcteq6x=(psbori(99,s,t,0, 0,11)+psbori(99,s,s-t,0, 0,11))*g1*g2

      return
      end

c-----------------------------------------------------------------------bg
      double precision function gglcteq6x(t,qt,xx1,xx2,s) !bg
c-----------------------------------------------------------------------------
c     gg->qq~
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s

      g1=cteq6(xx1,qt,0.,2,0)
      g2=cteq6(xx2,qt,0.,2,0)

      norm=min(3,naflav)
      gglcteq6x=(psbori(99,s,t,0, 0,10)+psbori(99,s,s-t,0, 0,10))
     .   *g1*g2/norm

      return
      end


c-----------------------------------------------------------------------bg
      double precision function cgcteq6x(t,qt,xx1,xx2,s) !bg
c-----------------------------------------------------------------------------
c     cg->cg
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s

      g1=cteq6(xx1,qt,0.,2,0)
      g2=cteq6(xx2,qt,0.,2,0)
      c1=cteq6(xx1,qt,0.,2,-4)       !c sea
      c2=cteq6(xx2,qt,0.,2,-4)       !c sea

      cgcteq6x=(psbori(99,s,t,0, 4,1)+psbori(99,s,s-t,0, 4,1))*
     *   (g2*c1+g1*c2)

      return
      end

c-----------------------------------------------------------------------bg
      double precision function qaqteq6x(t,qt,xx1,xx2,s) !bg
c-----------------------------------------------------------------------------
c     qq~->cc~
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s

      ub1=cteq6(xx1,qt,0.,2,-1)      !aq sea
      u1=cteq6(xx1,qt,0.,2,1)+ub1    !q val + q sea
      db1=cteq6(xx1,qt,0.,2,-2)      !aq sea
      d1=cteq6(xx1,qt,0.,2,2)+db1   !q val + q sea
      sb1=cteq6(xx1,qt,0.,2,-3)      !aq sea
      s1=sb1                         !q sea
      ub2=cteq6(xx2,qt,0.,2,-1)
      u2=cteq6(xx2,qt,0.,2,1)+ub2
      db2=cteq6(xx2,qt,0.,2,-2)
      d2=cteq6(xx2,qt,0.,2,2)+db2
      sb2=cteq6(xx2,qt,0.,2,-3)
      s2=sb2

      qaqteq6x=(psbori(99,s,t,0, 0,13)+psbori(99,s,s-t,0, 0,13))*
     *   (ub1*u2+ub2*u1+db1*d2+db2*d1+sb1*s2+sb2*s1)

      return
      end


c-----------------------------------------------------------------------bg
      double precision function qpcteq6x(t,qt,xx1,xx2,s) !bg
c-----------------------------------------------------------------------------
c     qc->qc q not charm or bottom
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s

      ub1=cteq6(xx1,qt,0.,2,-1)      !aq sea
      u1=cteq6(xx1,qt,0.,2,1)+ub1    !q val + q sea
      db1=cteq6(xx1,qt,0.,2,-2)      !aq sea
      d1=cteq6(xx1,qt,0.,2,2)+db1   !q val + q sea
      sb1=cteq6(xx1,qt,0.,2,-3)      !aq sea
      s1=sb1                         !q sea
      ub2=cteq6(xx2,qt,0.,2,-1)
      u2=cteq6(xx2,qt,0.,2,1)+ub2
      db2=cteq6(xx2,qt,0.,2,-2)
      d2=cteq6(xx2,qt,0.,2,2)+db2
      sb2=cteq6(xx2,qt,0.,2,-3)
      s2=sb2
      c1=cteq6(xx1,qt,0.,2,-4)       !c sea
      c2=cteq6(xx2,qt,0.,2,-4)       !c sea

      qpcteq6x=(psbori(99,s,t,1,4,1)+psbori(99,s,s-t,1,4,1))*
     *   ((ub1+u1+db1+d1+sb1+s1)*c2+
     *    (ub2+u2+db2+d2+sb2+s2)*c1)

      return
      end


csp**********************************************************************
csp                            GRV94
csp**********************************************************************
c------------------------------------------------------------------------
      function psdfh4(x,qqs,qq,icq,iq)
c------------------------------------------------------------------------
c psdfh4 - GRV structure functions (x*PDF)
c------------------------------------------------------------------------
      common /psar36/ alvc,qnorm,qnormp

      psdfh4=0.
      if(x.gt..99999)return
      if(icq.eq.2)then
        sq=log(log(qqs/.232**2)/log(.23/.232**2))
        if(iq.eq.0)then                                 !gluon
          alg=.524
          betg=1.088
          aag=1.742-.93*sq
          bbg=-.399*sq**2
          ag=7.486-2.185*sq
          bg=16.69-22.74*sq+5.779*sq*sq
          cg=-25.59+29.71*sq-7.296*sq*sq
          dg=2.792+2.215*sq+.422*sq*sq-.104*sq*sq*sq
          eg=.807+2.005*sq
          eeg=3.841+.361*sq
          psdfh4=(1.-x)**dg*(x**aag*(ag+bg*x+cg*x**2)*log(1./x)**bbg
     *    +sq**alg*exp(-eg+sqrt(eeg*sq**betg*log(1./x))))
        elseif(iq.eq.1.or.iq.eq.2)then !u_v or d_v
          if(sq.le.0d0)return  !single precision limit
          aau=.59-.024*sq
          bbu=.131+.063*sq
          auu=2.284+.802*sq+.055*sq*sq
          au=-.449-.138*sq-.076*sq*sq
          bu=.213+2.669*sq-.728*sq*sq
          cu=8.854-9.135*sq+1.979*sq*sq
          du=2.997+.753*sq-.076*sq*sq
          uv=auu*x**aau*(1.-x)**du*
     *    (1.+au*x**bbu+bu*x+cu*x**1.5)

          aad=.376
          bbd=.486+.062*sq
          add=.371+.083*sq+.039*sq*sq
          ad=-.509+3.31*sq-1.248*sq*sq
          bd=12.41-10.52*sq+2.267*sq*sq
          ccd=6.373-6.208*sq+1.418*sq*sq
          dd=3.691+.799*sq-.071*sq*sq
          dv=add*x**aad*(1.-x)**dd*
     *    (1.+ad*x**bbd+bd*x+ccd*x**1.5)

          aadel=.409-.005*sq
          bbdel=.799+.071*sq
          addel=.082+.014*sq+.008*sq*sq
          adel=-38.07+36.13*sq-.656*sq*sq
          bdel=90.31-74.15*sq+7.645*sq*sq
          ccdel=0.
          ddel=7.486+1.217*sq-.159*sq*sq
          delv=addel*x**aadel*(1.-x)**ddel*
     *    (1.+adel*x**bbdel+bdel*x+ccdel*x**1.5)

          alud=1.451
          betud=.271
          aaud=.41-.232*sq
          bbud=.534-.457*sq
          aud=.89-.14*sq
          bud=-.981
          cud=.32+.683*sq
          dud=4.752+1.164*sq+.286*sq*sq
          eud=4.119+1.713*sq
          eeud=.682+2.978*sq
          udsea=(1.-x)**dud*(x**aaud*(aud+bud*x+cud*x**2)
     *    *log(1./x)**bbud+sq**alud*exp(-eud+sqrt(eeud*sq**betud*
     *    log(1./x))))

          if(iq.eq.1)then                              !u_v
            psdfh4=uv
          elseif(iq.eq.2)then                          !d_v
            psdfh4=dv
          endif
        elseif(iq.eq.-3)then                           !s_sea
          als=.914
          bets=.577
          aas=1.798-.596*sq
          as=-5.548+3.669*sqrt(sq)-.616*sq
          bs=18.92-16.73*sqrt(sq)+5.168*sq
          ds=6.379-.35*sq+.142*sq*sq
          es=3.981+1.638*sq
          ees=6.402
          psdfh4=(1.-x)**ds*sq**als/log(1./x)**aas*(1.+as*sqrt(x)
     *    +bs*x)*exp(-es+sqrt(ees*sq**bets*log(1./x)))
        elseif(iabs(iq).lt.3)then                      !u_sea or d_sea
          aadel=.409-.005*sq
          bbdel=.799+.071*sq
          addel=.082+.014*sq+.008*sq*sq
          adel=-38.07+36.13*sq-.656*sq*sq
          bdel=90.31-74.15*sq+7.645*sq*sq
          ccdel=0.
          ddel=7.486+1.217*sq-.159*sq*sq
          delv=addel*x**aadel*(1.-x)**ddel*
     *    (1.+adel*x**bbdel+bdel*x+ccdel*x**1.5)

          alud=1.451
          betud=.271
          aaud=.41-.232*sq
          bbud=.534-.457*sq
          aud=.89-.14*sq
          bud=-.981
          cud=.32+.683*sq
          dud=4.752+1.164*sq+.286*sq*sq
          eud=4.119+1.713*sq
          eeud=.682+2.978*sq
          udsea=(1.-x)**dud*(x**aaud*(aud+bud*x+cud*x**2)
     *    *log(1./x)**bbud+sq**alud*exp(-eud+sqrt(eeud*sq**betud*
     *    log(1./x))))

          if(iq.eq.-1)then                           !u_sea
            psdfh4=(udsea-delv)/2.
          elseif(iq.eq.-2)then                       !d_sea
            psdfh4=(udsea+delv)/2.
          endif
        else
          psdfh4=0.
        endif
      elseif(icq.eq.1.or.icq.eq.3)then
        sq=log(log(qqs/.204**2)/log(.26/.204**2))
        if(iq.eq.1.or.iq.eq.2)then
          aapi=.517-.02*sq
          api=-.037-.578*sq
          bpi=.241+.251*sq
          dpi=.383+.624*sq
          anorm=1.212+.498*sq+.009*sq**2
          psdfh4=.5*anorm*x**aapi*(1.-x)**dpi*
     *    (1.+api*sqrt(x)+bpi*x)
        elseif(iq.eq.0)then
          alfpi=.504
          betpi=.226
          aapi=2.251-1.339*sqrt(sq)
          api=2.668-1.265*sq+.156*sq**2
          bbpi=0.
          bpi=-1.839+.386*sq
          cpi=-1.014+.92*sq-.101*sq**2
          dpi=-.077+1.466*sq
          epi=1.245+1.833*sq
          eppi=.51+3.844*sq
          psdfh4=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *    log(1./x)**bbpi+sq**alfpi*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
        elseif(iq.eq.-3)then
          alfpi=.823
          betpi=.65
          aapi=1.036-.709*sq
          api=-1.245+.713*sq
          bpi=5.58-1.281*sq
          dpi=2.746-.191*sq
          epi=5.101+1.294*sq
          eppi=4.854-.437*sq
          psdfh4=sq**alfpi/log(1./x)**aapi*(1.-x)**dpi*
     *    (1.+api*sqrt(x)+bpi*x)*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x)))
        elseif(iabs(iq).lt.3)then
          alfpi=1.147
          betpi=1.241
          aapi=.309-.134*sqrt(sq)
          api=.219-.054*sq
          bbpi=.893-.264*sqrt(sq)
          bpi=-.593+.24*sq
          cpi=1.1-.452*sq
          dpi=3.526+.491*sq
          epi=4.521+1.583*sq
          eppi=3.102
          psdfh4=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *    log(1./x)**bbpi+sq**alfpi*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
        else
          psdfh4=0.
        endif
      elseif(icq.eq.0)then
        sq=log(log(qqs/.204**2)/log(.26/.204**2))
        if(iq.eq.0)then
          alfpi=.504
          betpi=.226
          aapi=2.251-1.339*sqrt(sq)
          api=2.668-1.265*sq+.156*sq**2
          bbpi=0.
          bpi=-1.839+.386*sq
          cpi=-1.014+.92*sq-.101*sq**2
          dpi=-.077+1.466*sq
          epi=1.245+1.833*sq
          eppi=.51+3.844*sq
          psdfh4=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *    log(1./x)**bbpi+sq**alfpi*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
     *    *.543
        else
          alfpi=.823
          betpi=.65
          aapi=1.036-.709*sq
          api=-1.245+.713*sq
          bpi=5.58-1.281*sq
          dpi=2.746-.191*sq
          epi=5.101+1.294*sq
          eppi=4.854-.437*sq
          str=sq**alfpi/log(1./x)**aapi*(1.-x)**dpi*
     *    (1.+api*sqrt(x)+bpi*x)*
     *    exp(-epi+sqrt(eppi*sq**betpi*log(1./x)))
          if(iq.eq.3)then
            psdfh4=str*.543*2.
          else
            aapi=.517-.02*sq
            api=-.037-.578*sq
            bpi=.241+.251*sq
            dpi=.383+.624*sq
            anorm=1.212+.498*sq+.009*sq**2
            val=.5*anorm*x**aapi*(1.-x)**dpi*
     *      (1.+api*sqrt(x)+bpi*x)

            alfpi=1.147
            betpi=1.241
            aapi=.309-.134*sqrt(sq)
            api=.219-.054*sq
            bbpi=.893-.264*sqrt(sq)
            bpi=-.593+.24*sq
            cpi=1.1-.452*sq
            dpi=3.526+.491*sq
            epi=4.521+1.583*sq
            eppi=3.102
            sea=(1.-x)**dpi*(x**aapi*(api+bpi*sqrt(x)+cpi*x)*
     *      log(1./x)**bbpi+sq**alfpi*
     *      exp(-epi+sqrt(eppi*sq**betpi*log(1./x))))
            if(iq.eq.1)then
              psdfh4=(.836*(val+2.*sea)-.587*str)
            elseif(iq.eq.2)then
              psdfh4=(.25*(val+2.*sea)+.587*str)
            else
              psdfh4=0.
            endif
          endif
        endif
        psdfh4=psdfh4/(1.+qq/.59)**2

      elseif(icq.eq.4.and.iq.eq.1)then
        psdfh4=x**3*(1.-x)**alvc*(alvc+3.)*(alvc+2.)*(alvc+1.)
      else
        psdfh4=0.
      endif
      return
      end

c-----------------------------------------------------------------------
      double precision function psjvrx(t,qt,xx1,xx2,s)
c-----------------------------------------------------------------------
#include "sem.h"
      double precision psbori,t,s,g1,ub1,u1,db1,d1,sb1,s1
     *                           ,g2,ub2,u2,db2,d2,sb2,s2

      g1=dble(psdfh4(xx1,qt,0.,2,0))
      ub1=dble(psdfh4(xx1,qt,0.,2,-1))
      u1=dble(psdfh4(xx1,qt,0.,2,1))+ub1
      db1=dble(psdfh4(xx1,qt,0.,2,-2))
      d1=dble(psdfh4(xx1,qt,0.,2,2))+db1
      sb1=dble(psdfh4(xx1,qt,0.,2,-3))
      s1=sb1
      g2=dble(psdfh4(xx2,qt,0.,2,0))
      ub2=dble(psdfh4(xx2,qt,0.,2,-1))
      u2=dble(psdfh4(xx2,qt,0.,2,1))+ub2
      db2=dble(psdfh4(xx2,qt,0.,2,-2))
      d2=dble(psdfh4(xx2,qt,0.,2,2))+db2
      sb2=dble(psdfh4(xx2,qt,0.,2,-3))
      s2=sb2

      psjvrx=g1*g2*(psbori(99,s,t,0,0,1)+psbori(99,s,s-t,0,0,1)
     *+psbori(99,s,t,0,0,2)+psbori(99,s,s-t,0,0,2))/2. !bg psbori changed. factor /2. included

     *+(psbori(99,s,t,0,1,1)+psbori(99,s,s-t,0,1,1))*
     *(g2*(u1+ub1+d1+db1+s1+sb1)+g1*(u2+ub2+d2+db2+s2+sb2))

     *+(psbori(99,s,t,1,1,1)+psbori(99,s,s-t,1,1,1))/2.* !bg
     *(u1*u2+ub1*ub2+d1*d2+db1*db2+s1*s2+sb1*sb2)

     *+(psbori(99,s,t,1,-1,1)+psbori(99,s,s-t,1,-1,1)
     *+psbori(99,s,t,1,-1,2)+
     *psbori(99,s,s-t,1,-1,2)+psbori(99,s,t,1,-1,3)
     *+psbori(99,s,s-t,1,-1,3))*
     *(u1*ub2+ub1*u2+d1*db2+db1*d2+s1*sb2+sb1*s2)

     *+(psbori(99,s,t,1,2,1)+psbori(99,s,s-t,1,2,1))*
     *((u1+ub1)*(d2+db2+s2+sb2)+(u2+ub2)*(d1+db1+s1+sb1)+
     *(d1+db1)*(u2+ub2+s2+sb2)+(d2+db2)*(u1+ub1+s1+sb1)+
     *(s1+sb1)*(u2+ub2+d2+db2)+(s2+sb2)*(u1+ub1+d1+db1))

      return
      end


csp***************************************************************************
csp           DUKE OWENS (very old pdf)
csp****************************************************************************

c-----------------------------------------------------------------------
      double precision function psjwox(t,qt,xx1,xx2,s)
c-----------------------------------------------------------------------
      double precision x,scale,upv1,dnv1,sea1,str1,chm1,gl1,
     *upv2,dnv2,sea2,str2,chm2,gl2,psbori,s,t
      scale=sqrt(qt)
      x=xx1
      call strdo1(x,scale,upv1,dnv1,sea1,str1,chm1,gl1)
      x=xx2
      call strdo1(x,scale,upv2,dnv2,sea2,str2,chm2,gl2)

      psjwox=gl1*gl2*(psbori(99,s,t,0,0,1)+psbori(99,s,s-t,0,0,1)
     *+psbori(99,s,t,0,0,2)+psbori(99,s,s-t,0,0,2)+psbori(99,s,t,0,0,3)
     *+psbori(99,s,s-t,0,0,3))/2. !bg
     *+(psbori(99,s,t,0,1,1)+psbori(99,s,s-t,0,1,1)
     *+psbori(99,s,t,0,1,2)+psbori(99,s,s-t,0,1,2)+psbori(99,s,t,0,1,3)
     *+psbori(99,s,s-t,0,1,3))*(gl2*(upv1+dnv1+4.*sea1+2.*str1+2.*chm1)+
     *gl1*(upv2+dnv2+4.*sea2+2.*str2+2.*chm2))
     *+(psbori(99,s,t,1,1,1)+psbori(99,s,s-t,1,1,1)
     *+psbori(99,s,t,1,1,2)+psbori(99,s,s-t,1,1,2)+psbori(99,s,t,1,1,3)+
     *psbori(99,s,s-t,1,1,3))/2.* !bg
     *((upv1+sea1)*(upv2+sea2)+(dnv1+sea1)*(dnv2+sea2)+2.*sea1*sea2
     *+2.*str1*str2+2.*chm1*chm2)
     *+(psbori(99,s,t,1,-1,1)+psbori(99,s,s-t,1,-1,1)
     *+psbori(99,s,t,1,-1,2)+
     *psbori(99,s,s-t,1,-1,2)+psbori(99,s,t,1,-1,3)
     *+psbori(99,s,s-t,1,-1,3))*
     *((upv1+sea1)*sea2+sea1*(upv2+sea2)+(dnv1+sea1)*sea2+
     *sea1*(dnv2+sea2)+2.*str1*str2+2.*chm1*chm2)
     *+(psbori(99,s,t,1,2,1)
     *+psbori(99,s,s-t,1,2,1)+psbori(99,s,t,1,2,2)
     *+psbori(99,s,s-t,1,2,2)
     *+psbori(99,s,t,1,2,3)+psbori(99,s,s-t,1,2,3))*
     *(upv1*dnv2+upv2*dnv1+(upv1+dnv1)*(2.*sea2+2.*str2+2.*chm2)+
     *(upv2+dnv2)*(2.*sea1+2.*str1+2.*chm1)+
     *4.*sea1*(2.*sea2+2.*str2+2.*chm2)+2.*str1*(4.*sea2+2.*chm2)+
     *2.*chm1*(4.*sea2+2.*str2))
      return
      end
c------------------------------------------------------------------------
      subroutine strdo1(x,scale,upv,dnv,sea,str,chm,gl)
c------------------------------------------------------------------------
c :::::::::::: duke owens set 1 ::::::::::::::::::::::::::::
c------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      double precision
     +       f(5),a(6,5),b1(3,6,5)
      data q0,ql1/2.d0,.2d0/
      data b1/3.d0,0.d0,0.d0,.419d0,.004383d0,-.007412d0,
     &3.46d0,.72432d0,-.065998d0,4.4d0,-4.8644d0,1.3274d0,
     &6*0.d0,1.d0,
     &0.d0,0.d0,.763d0,-.23696d0,.025836d0,4.d0,.62664d0,-.019163d0,
     &0.d0,-.42068d0,.032809d0,6*0.d0,1.265d0,-1.1323d0,.29268d0,
     &0.d0,-.37162d0,-.028977d0,8.05d0,1.5877d0,-.15291d0,
     &0.d0,6.3059d0,-.27342d0,0.d0,-10.543d0,-3.1674d0,
     &0.d0,14.698d0,9.798d0,0.d0,.13479d0,-.074693d0,
     &-.0355d0,-.22237d0,-.057685d0,6.3494d0,3.2649d0,-.90945d0,
     &0.d0,-3.0331d0,1.5042d0,0.d0,17.431d0,-11.255d0,
     &0.d0,-17.861d0,15.571d0,1.564d0,-1.7112d0,.63751d0,
     &0.d0,-.94892d0,.32505d0,6.d0,1.4345d0,-1.0485d0,
     &9.d0,-7.1858d0,.25494d0,0.d0,-16.457d0,10.947d0,
     &0.d0,15.261d0,-10.085d0/
      wn=1.d0
      s= log( log( max(q0,scale)/ql1)/ log(q0/ql1))
      do 10 i=1,5
      do 10 j=1,6
   10 a(j,i)=b1(1,j,i)+s*(b1(2,j,i)+s*b1(3,j,i))
      do 40 i=1,5
   40 f(i)=a(1,i)*x**a(2,i)*(wn-x)**a(3,i)*(wn+x*
     &    (a(4,i)+x*(a(5,i)+x*a(6,i))))
      do 50 i=1,2
      aa=wn+a(2,i)+a(3,i)
   50 f(i)=f(i)*utgam2(aa)/((wn+a(2,i)*a(4,i)/aa)
     &*utgam2(a(2,i))*utgam2(wn+a(3,i)))
      upv=f(1)-f(2)
      dnv=f(2)
      sea=f(3)/6.d0
      str=sea
      chm=f(4)
      gl =f(5)
      return
      end

c------------------------------------------------------------------------
      function fractalpdf(x,qqs,qq,icq,iq)
c------------------------------------------------------------------------
c fractalpdf - structure functions based on fit by Lastovicka (DESY 02-036,
c              March 2002)) fit. Fit function based on fractal dimension of
c              proton structure function.
c Global fit (qq, iq and icq are dummy)
c------------------------------------------------------------------------

c self-similarity of proton structure function only for x<0.01 (no valence)
      if(x.gt.0.9)then
        fractalpdf=0.
        return
      endif

      iq0=iq
      iqc0=icq
      qq0=qq
c fit parameters:
      D0=0.523
      D1=0.074
      D2=1.
      D3=-1.282
      Q02=0.051
c scaling variable
      Qmag=1.+qqs/Q02

c function
      fractalpdf=exp(D0)*Q02*x**(-D2+1.)/(1.+D3-D1*log(x))
     &          *(x**(-D1*log(Qmag))*Qmag**(D3+1.)-1.)

      return
      end

