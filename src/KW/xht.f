C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

  
  
!  message for Benjamin :
!
!   z z s o f t   is not defined used more,
!   
!   so I replace it by     0.1199000  
!   
!        look for  !kw
   
   


c------------------------------------------------------------------------!bg
      subroutine testpsahot(k,n)                               !bg
c------------------------------------------------------------------------!bg
#include "aaa.h"
#include "ems.h"
#include "sem.h"
      integer cpt(110) !counter
      character mod1*1, mod2*1,mod3*1,mod4*7,mod5*5
      common/xpit/xp1t,xp2t
      common/iipomend/ipomend,iqqt,cpt

      ipomend=1    !activate test inside psahot
      iret=0
      delta=0.01

      write(ifmt,'(72a)')'+',('-',i=1,70),'+'
      write(ifmt,'(a,30x,a,30x,a)')'|','testpashot','|'
      write(ifmt,'(72a)')'+',('-',i=1,70),'+'

       write(ifhi,'(a)')       'orientation landscape'
       write(ifhi,'(a)')       'newpage'
       write(ifhi,'(a)')'zone 3 3 1 set scalel 0.5  set scalem 0.6'
       write(ifhi,'(a)')       'set basl 0.004'
       write(ifhi,'(a)')        'cd .'

      do itypep=1,4
      iqqt=itypep-1
      write(ifmt,'(1x,a,i1,a,$)')'iqq',iqqt,' '
      xp1t=0.1
      do while(xp1t.le.0.9)
        xp2t=0.1
        do while(xp2t.le.0.9)
          do i=1,110
            cpt(i)=0
          enddo
ctp20170105 value of q2min given in ems.f in q2cmin
          do i=1,10000
            call psahot(k,n,iret)
          enddo
          norm=0 !for normalisation
          do i=1,110
          norm=norm+cpt(i)
          enddo

          kmod=aint(xp1t*10.)
          lmod=aint(xp2t*10.)
          write(mod1,'(I1)')kmod
          write(mod2,'(I1)')lmod
          write(mod3,'(I1)')iqqt
          mod5=mod3//mod1//mod2
          if(iqqt.eq.0) then
            mod4='sea-sea'
          elseif(iqqt.eq.1) then
            mod4='val-sea'
          elseif(iqqt.eq.2) then
            mod4='sea-val'
          else
            mod4='val-val'
          endif

          write(ifhi,'(a)')       '!----------------------------------'
          write(ifhi,'(a)')       '!   Pomeron xp distribution '
          write(ifhi,'(a)')       '!----------------------------------'

          write(ifhi,'(a)')  'openhisto name xPomMC'//mod5
          write(ifhi,'(a)')  'htyp lru xmod lin ymod log'
          write(ifhi,'(a)')  'xrange 0. 1. yrange auto auto'
          write(ifhi,'(a)')  'text 0 0 "xaxis xp"'
          write(ifhi,'(a)')  'text 0 0 "yaxis F'//mod4//'"'
          write(ifhi,'(a)')  'text   0.5   0.7  "xp1='//mod1//
     *                       'xp2='//mod2//'" '
          write(ifhi,'(a)')  'histoweight 1'
          write(ifhi,'(a)')  ' array 2'
          do i=1,110
            write(ifhi,*) (i-1)*delta,real(cpt(i))/real(norm)/delta
          enddo
          write(ifhi,'(a)')'endarray closehisto plot 0-'
          call distributionX(xp1t,xp2t,iqqt,n,k)
          enddo
        xp2t=xp2t+0.2
        enddo
        xp1t=xp1t+0.2
        enddo
        ipomend=0
        write(ifmt,'(a)')' testpsahot finished'
        stop '\n\n    STOP at the end of testpsahot\n\n'
        end

c-----------------------------------------------------------------------!bg
      subroutine distributionX(xp1,xp2,iqq,ncolp,kcol)              !bg
c-----------------------------------------------------------------------!bg
c
c     theoretical distribution of xp for given pomeron ends xp1,xp2
c     xp is all cases for the valence side
c-----------------------------------------------------------------------
      common /psar7/  delx,alam3p,gam3p
#include "aaa.h"
#include "sem.h"
#include "ems.h"
      double precision xp,xm,EsaturGluonTil,EsaturQuarkTil
     *                ,EsaturValTil
      character mod1*5, mod2*5,mod3*5
      dimension qq(2)

      qq(1)=q2cmin(1)
      qq(2)=q2cmin(2)

       kmod=aint(xp1*10)
       lmod=aint(xp2*10)
       write(mod1,'(I1)')kmod
       write(mod2,'(I1)')lmod
       write(mod3,'(I1)')iqq
       mod3=mod3(1:1)//mod1(1:1)//mod2(1:1)

       write(ifhi,'(a)')  'openhisto name xPomth'//mod3
       write(ifhi,'(a)')  'htyp lru xmod lin ymod log'
       write(ifhi,'(a)')  'xrange 0. 1. yrange auto auto'
       write(ifhi,'(a)')  'text 0 0 "xaxis xp"'
       write(ifhi,'(a)')  'text 0 0 "yaxis F"'
       write(ifhi,'(a)')  'histoweight 1'
       write(ifhi,'(a)')  ' array 2'
c-------definition of parameters-----c
      !b_nk = bhpr(ncolp,kcol) !presently comented in ems.h
      ndummy= ncolp
      b_k = bk(kcol)
      b = b_k   !b_nk
      sy=xp1*xp2*engy*engy
      je1=2
      je2=2 !all emissions
      delta=0.01
      spmin=4.*max(q2cmin(1),q2cmin(2)) !??????????????????????
      if(sy.le.spmin)goto 998
      fnorm=om52pp(sy,xp1,xp2,b,iqq,je1,je2) ! for normalisation !tp10142013:factor 2 or 0.5 ???
c---------------------------------c

      ef1=0
      ef2=0
      ef3=0
      ef4=0
      if( je1.ge.1             .and. je2.ge.1)             ef1=1
      if( je1.ge.1             .and.(je2.eq.0.or.je2.eq.2))ef2=1
      if((je1.eq.0.or.je1.eq.2).and. je2.ge.1)             ef3=1
      if((je1.eq.0.or.je1.eq.2).and.(je2.eq.0.or.je2.eq.2))ef4=1


      iclv=2   !PDF fitted on proton, so deconvolution done for nucleons
      if(iqq.eq.2)then
        icl=icltar
      else
        icl=iclpro
      endif
      xppii=0d0
      if(iqq.eq.1)then
        xppii=xp1
      elseif(iqq.eq.2)then
        xppii=xp2
      endif

      delss=0.5*(dels(1)+dels(2))               !just to optimize integration
      if(iqq.eq.3)delss=-0.5   !better value for hard (no soft evolution)

      xmin=spmin/sy

      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
      ximin=0.
      !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

c      alpq=0.5*(betff(1)+betff(2))
      xma1=1.

      !----------------------------------------------------------
      !the following is used used to compute (in the loops)
      !the factor g(b) for the b dependence, without 1/.0389:
      !g'(b)= 1/(4*pi*r2hhx) * exp(-b**2/(4.*.0389*r2hhx))
      !where r2hhx depends on zh or xm .
      !----------------------------------------------------------
      r2hh=r2had(iclpro)+r2had(icltar)
      r2hhs=r2hh+slopom*log(max(1.,sy))
      zbb=exp(-b**2/(4.*.0389*r2hhs))
      fopi=4.*pi

c--------bounds for xp--------------------
      uumin=0.0001
      uumax=1.
      if(iqq.eq.0)uumin=0.1199000  !kw
c----------------------------------------
      uu=uumin
      do while(uu.le.uumax)
        res=0.
        inorm=0
c--------bounds for xm--------------------
      vvmin=xmin/uu
      if(vvmin.ge.1.)goto 1
      vvmax=1.
      if(iqq.ne.3)then
        if(vvmin.lt.0.1199000)vvmin=0.1199000  !kw
      endif
c----------------------------------------
      vv=vvmin
      do while(vv.le.vvmax)
        stq=0.
        zh=uu*vv
        xp=dble(uu)
        xm=dble(vv)
        if(zh.gt.xmin) then
          sgg=       ef1  *pijet(2,zh*sy,0,0)
     *         + (ef2+ef3)*pijet(1,zh*sy,0,0)
     *         +     ef4  *pijet(0,zh*sy,0,0)
          sgq=       ef1  *pijet(2,zh*sy,0,1)
     *         +     ef2  *pijet(1,zh*sy,1,0)
     *         +     ef3  *pijet(1,zh*sy,0,1)
     *         +     ef4  *pijet(0,zh*sy,0,1)
          sqq=       ef1  *pijet(2,zh*sy,1,1)
     *         + (ef2+ef3)*pijet(1,zh*sy,1,1)
     *         +     ef4  *pijet(0,zh*sy,1,1)
          sqaq=      ef1  *pijet(2,zh*sy,-1,1)
     *         + (ef2+ef3)*pijet(1,zh*sy,-1,1)
     *         +     ef4  *pijet(0,zh*sy,-1,1)
          sqqp=      ef1  *pijet(2,zh*sy,1,2)
     *         + (ef2+ef3)*pijet(1,zh*sy,1,2)
     *         +     ef4  *pijet(0,zh*sy,1,2)

          sqqi=sqq

          if(iqq.eq.0) then
                glu1=sngl(EsaturGluonTil(xp,qq(1),1))
                sea1=sngl(EsaturQuarkTil(xp,qq(1),1,999)) 
     .                     *2*noflav(qq(1))
                glu2=sngl(EsaturGluonTil(xm,qq(2),2))
                sea2=sngl(EsaturQuarkTil(xm,qq(2),2,999)) 
     .                     *2*noflav(qq(2)) 
                r2hhx=r2hh-slopom*log(zh)
                stq= glu1*glu2*sgg
     *               +(glu1*sea2+sea1*glu2)*sgq
     *               +sea1*sea2*sqq
                stq=stq * zbb**(r2hhs/r2hhx) /r2hhx/fopi  !g(b)
                stq=stq*zh**(-1-delss)    !this factor cancel with the jacobian in om52pp
          endif

          if(iqq.eq.1.or.iqq.eq.2) then
            if(xp*dble(xppii).lt..99999d0)then
              uv1=sngl(EsaturValTil(1,xp,dble(xppii),qq(iqq),1,icl))
              dv1=sngl(EsaturValTil(1,xp,dble(xppii),qq(iqq),2,icl))
              glu2=sngl(EsaturGluonTil(xm,qq(3-iqq),3-iqq))
              sea2=sngl(EsaturQuarkTil(xm,qq(3-iqq),3-iqq,999))  
     .                  *2*noflav(qq(3-iqq))
              r2hhx=r2hh-slopom*sngl(log(xm))
              if(xp.ne.1.)
     *             stq=(glu2*sgq+sea2*sqq)*(uv1+dv1)*zh**(-1-delss)
              stq=stq * zbb**(r2hhs/r2hhx) /r2hhx/fopi !g(b)
            endif
          endif

          if(iqq.eq.3) then
           if(xp*dble(xp1).le..9999d0.and.xm*dble(xp2).le..9999d0
     *     .or.xm*dble(xp1).le..9999d0
     *     .and.xp*dble(xp2).le..9999d0)then
            stq=
     *      (psharf(xp,dble(xp1),xm,dble(xp2),sqqi,sqqp,sqaq))
             stq=stq*zbb**(r2hhs/r2hh)/r2hh/fopi  !dep in b.
           endif
          endif

          stq=stq* factk / sigine * 10 / 2.
          res=res+stq*delta !for integration over xm
          inorm=inorm+1
        endif
      vv=vv+delta
      enddo
        if(inorm.ne.0) then
        write(ifhi,*) sngl(xp),res/fnorm
        endif
    1 continue
      uu=uu+delta
      enddo

      write(ifhi,'(a)')'endarray closehisto plot 0'
 998  continue
      end

