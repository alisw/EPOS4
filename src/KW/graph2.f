C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c#######################################################################
c#######################################################################
c################## Output data as histo file ##########################
c#######################################################################
c#######################################################################

c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> create a text file from the input options file name with the suffix 
! histo2
!
!-----------------------------------------------------------------------
      subroutine histo2BeginFigure()
c#######################################################################
#include "aaa.h"
      character :: cc*5, txt*80
      integer :: i, ii, iii, ii1, icanvas, m, n
      data icanvas /0/
      save icanvas

      icanvas=icanvas+1

      ii=0
      if(icanvas.eq.1)then
         nopenr=nopenr+1
         ii=index(fnhi(1:nfnhi),".histo")-1

         write(cc,'(i5)')icanvas
         n=5-(int(log10(float(icanvas)))+1)
         do m=2,n
            cc(m:m)='0'
         enddo
         cc(1:1)='_'

         ii1=0
         do i=1,ii
            if(fnhi(i:i).eq.'/') ii1=i
         enddo

         ii1=ii1+1
         iii=ii-ii1+1
         txt(1:iii)=fnhi(ii1:ii)
         do m=1,iii
            if(txt(m:m).eq.'-')txt(m:m)='_'
         enddo

         open(97,file=fnhi(1:ii)//".histo2")

      endif

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write zone informations
!
!> @param[in] nb_rows number of rows
!> @param[in] nb_cols number of columns
!-----------------------------------------------------------------------
      subroutine histo2BeginPlot(nb_rows, nb_cols)
c#######################################################################
      integer :: nb_rows, nb_cols

      write(97,'(a, i4, i4, a)') 
     .        'zone', nb_rows, nb_cols, ' 1'

      end
      

c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write histo informations
!
!> @param[in] name name of the calling subroutine
!> @param[in] text test of the histo
!-----------------------------------------------------------------------
      subroutine histo2BeginSubPlot(name, text)
c#######################################################################
      character(*) :: name, text

      write(97,'(3a/3a)')
     .     'openhisto2 name "', text, '"',
     .     'subroutine "', name, '"'

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write histo informations
!
!> @param[in] title test of the histo
!> @param[in] xmin min value in x
!> @param[in] xmax max value in x
!> @param[in] xstep step in x
!> @param[in] ymin min value in y
!> @param[in] ymax max value in y
!> @param[in] ystep step in y
!> @param[in] xtitle title for x-axis
!> @param[in] ytitle title for y-axis
!-----------------------------------------------------------------------
      subroutine histo2AddHeader(title,
     .     xmin, xmax, xstep,
     .     ymin, ymax, ystep,
     .     xtitle, ytitle)
c#######################################################################
      character(*) :: title, xtitle, ytitle
      real :: xmin, xmax, xstep, ymin, ymax, ystep

      write(97,'(3a/(a, 2e11.3)/(a, 2e11.3)/(a, 2e11.3)/3a/3a)')
     .     'title "', title, '"',
     .     'xrange ', xmin, xmax, 
     .     'yrange ', ymin, ymax,
     .     'step ', xstep, ystep,
     .     'txt "xaxis ', xtitle, '"',
     .     'txt "yaxis ', ytitle, '"'
      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write histo informations related to PHSD
!
!-----------------------------------------------------------------------
      subroutine histo2AddPHSDHeader()
c#######################################################################
      integer iphsd
      common/cphsd/iphsd
      real    egymin,egymax,elab,ecms,ekin
      common/enrgy/egymin,egymax,elab,ecms,ekin
      real         bmaxim,bminim,phimax,phimin,zzsoft,asatur,bsatur
      common/nucl2/bmaxim,bminim,phimax,phimin,zzsoft,asatur,bsatur
      integer      nevent,nfull,nfreeze,ninicon
      common/events/nevent,nfull,nfreeze,ninicon
      character(len = 70) :: text1 
      character(len = 15) :: text2, text3
      character(len = 10) :: energy
      character(len = 25) :: bstring

      if (iphsd.eq.2) then
         write(text2,'(a,i6)')'num=',nfreeze
         text3 = "EPOSi + PHSDe"
      elseif (iphsd.eq.9) then
         write(text2,'(a,i6)')'num=',ninicon
         text3 = "pure PHSD"
      endif

      write(energy,'(f8.2)')ecms
      
      if (bminim.eq.bmaxim) then
         write(bstring,'(a,f5.2,a)')'b=',bminim,'fm'
      else
         write(bstring,'(f5.2,a,f5.2,a)')bminim,'fm<=b<=',bmaxim,'fm'
      endif

      text1 = "AuAu@"//energy//"GeV, "//bstring//", "//text2
      call stripSpaces(text1)

      write(97,'(3a/3a/3a)')
     .     'text  0.05 0.93  " ', trim(text1), ' "',
     .     'text  0.05 0.05  " ', trim(text3), ' "'
      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write array data
!
!> @param[in] ndim dimension
!> @param[in] x x value
!> @param[in] y y value
!> @param[in] z z value
!-----------------------------------------------------------------------
      subroutine histo2FillArray(ndim, x, y, z)
c#######################################################################
      integer :: ndim
      real :: x, y, z

      if((ndim.eq.2) .and. (z.ne.0)) then
         write(97,'(2(e11.3, a), e11.3)')
!-----------------------------------------------------------------
     .        x, ' ', y, ' ', z
!-----------------------------------------------------------------
      endif

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write array data
!
!> @param[in] x_center
!> @param[in] y_center
!> @param[in] radius
!-----------------------------------------------------------------------
      subroutine histo2DrawCircle(x_center, y_center, radius)
      real :: x_center, y_center, radius
      
      write(97,'(3(a, e11.3))') 
     .     'circle ', x_center, ' ', y_center, ' ', radius

      end



c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write array footer
!
!-----------------------------------------------------------------------
      subroutine histo2BeginArray()
c#######################################################################

      write(97,'(a)') 'array 3'

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write array footer
!
!-----------------------------------------------------------------------
      subroutine histo2EndArray()
c#######################################################################

      write(97,'(a)') 'endarray'

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write histo footer
!
!-----------------------------------------------------------------------
      subroutine histo2EndSubPlot()
c#######################################################################

      write(97,'(a)') 'closehisto2'

      end


c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write plot command
!
!-----------------------------------------------------------------------
      subroutine histo2EndPlot()
c#######################################################################
      write(97,'(a)') 'plot 0-'

      end

      
c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> write plot command
!
!-----------------------------------------------------------------------
      subroutine histo2EndLastPlot()
c#######################################################################
      write(97,'(a)') 'plot 0'

      end

      
c#######################################################################
!-----------------------------------------------------------------------
!> @brief
!> close histo file
!
!-----------------------------------------------------------------------
      subroutine histo2EndFigure()
c#######################################################################
#include "aaa.h"

      if(nopenr.gt.0)then
         close (97)
      endif

      end



