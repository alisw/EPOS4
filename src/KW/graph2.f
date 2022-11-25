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
      subroutine histo2BeginSubPlot(name, text, title,
     .     xmin, xmax, xstep,
     .     ymin, ymax, ystep,
     .     xtitle, ytitle)
c#######################################################################
      character(*) :: name, text, title, xtitle, ytitle
      real :: xmin, xmax, xstep, ymin, ymax, ystep

      write(97,'(3a/3a/3a/(a, 2e11.3)/(a, 2e11.3)/(a, 2e11.3)
     .           /3a/3a/a)')
     .     'openhisto2 name "', text, '"',
     .     'subroutine "', name, '"',
     .     'title "', title, '"',
     .     'xrange ', xmin, xmax, 
     .     'yrange ', ymin, ymax,
     .     'step ', xstep, ystep,
     .     'txt "xaxis ', xtitle, '"',
     .     'txt "yaxis ', ytitle, '"',
     .     'array 3'

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
!> write histo footer
!
!-----------------------------------------------------------------------
      subroutine histo2EndSubPlot()
c#######################################################################

      write(97,'(a)') 'endarray'
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
      write(97,'(a)') 'plot'

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



