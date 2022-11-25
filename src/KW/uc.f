C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

      subroutine readcharm
      implicit none
      include 'coms.f'

      character*500  fnin,fnus1,fnus2,fnus3,fnus4
      common/fname2/ fnin,fnus1,fnus2,fnus3,fnus4
      
      real*8 dummy
      integer icm
      

c     read in charm cross sections-------------
      open(UNIT=81,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Dpi_cross_secs_el.dat')
      do icm = 1,1000
         read(81,*)dummy,Dppip(icm),Dppi0(icm),dummy,Dppim(icm),dummy
     $        ,dummy,Dzpip(icm),dummy,dummy,Dzpim(icm)
      enddo
      close(UNIT=81)

      open(UNIT=82,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Dstarpi_cross_secs_el.dat')
      do icm = 1,1000
         read(82,*)dummy,Dsppip(icm),Dsppi0(icm),dummy,Dsppim(icm),dummy
     $        ,dummy,Dszpip(icm),dummy,dummy,Dszpim(icm)
      enddo
      close(UNIT=82)

      open(UNIT=83,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'DK_cross_secs_el.dat')
      do icm = 1,1000
         read(83,*)dummy,DpKp(icm),DpKz(icm),dummy,dummy
     $        ,DzKp(icm),DzKz(icm)
      enddo
      close(UNIT=83)

      open(UNIT=84,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'DKbar_cross_secs_el.dat')
      do icm = 1,1000
         read(84,*)dummy,DpaKz(icm),DpKm(icm),dummy,dummy
     $        ,DzaKz(icm),DzKm(icm)
      enddo
      close(UNIT=84)

      open(UNIT=85,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'DstarK_cross_secs_el.dat')
      do icm = 1,1000
         read(85,*)dummy,DspKp(icm),DspKz(icm),dummy,dummy
     $        ,DszKp(icm),DszKz(icm)
      enddo
      close(UNIT=85)

      open(UNIT=86,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'DstarKbar_cross_secs_el.dat')
      do icm = 1,1000
         read(86,*)dummy,DspaKz(icm),DspKm(icm),dummy,dummy
     $        ,DszaKz(icm),DszKm(icm)
      enddo
      close(UNIT=86)

      open(UNIT=87,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Deta_cross_secs_el.dat')
      do icm = 1,1000
         read(87,*)dummy,dummy,Deta(icm)
      enddo
      close(UNIT=87)

      open(UNIT=88,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Dstareta_cross_secs_el.dat')
      do icm = 1,1000
         read(88,*)dummy,dummy,Dseta(icm)
      enddo
      close(UNIT=88)

      open(UNIT=89,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'dn_el.dat')
      do icm = 1,7901
         read(89,*)dummy,Dzn(icm),dummy,Dzp(icm),dummy,dummy,dummy
      enddo
      close(UNIT=89)

      open(UNIT=90,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'dbarn_el.dat')
      do icm = 1,7909
         read(90,*)dummy,Dmn(icm),dummy,Dmp(icm),dummy,dummy,dummy
      enddo
      close(UNIT=90)

      open(UNIT=91,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Ddelta_cross_secs_el.dat')
      do icm = 1,5000
         read(91,*)dummy,DpDelpp(icm),DpDelp(icm),dummy,DpDelz(icm),
     $        dummy,DpDelm(icm),dummy,dummy,dummy,dummy,dummy,dummy
     $        ,dummy,dummy
      enddo
      close(UNIT=91)

      open(UNIT=92,
     $FILE=fnus3(1:index(fnus3,' ')-1)//'Dbardelta_cross_secs_el.dat')
      do icm = 1,5204
         read(92,*)dummy,DbzDelpp(icm),DbzDelp(icm),dummy,DbzDelz(icm),
     $        dummy,DbzDelm(icm),dummy,dummy,dummy,dummy,dummy,dummy
     $        ,dummy,dummy
      enddo
      close(UNIT=92)


      return
      end
