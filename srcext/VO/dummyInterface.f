C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine hqini
!-----------------------------------------------------------------------
      subroutine hqini
      end subroutine hqini


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine exiteposevent
!
!> @param[in] i
!> @param[in] j
!> @param[in] k
!-----------------------------------------------------------------------
      subroutine exiteposevent(i,j,k)
      integer i, j, k
      end subroutine exiteposevent


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine aacharm
!-----------------------------------------------------------------------
      subroutine aacharm
      end subroutine aacharm


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine hqoptns
!
!> @param[in] line
!> @param[in] j
!-----------------------------------------------------------------------
      subroutine hqoptns(line,j)
      character*1000 line
      integer j
      end subroutine hqoptns


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for function idKlaus
!
!> @param[in] idPB
!> @return 0
!-----------------------------------------------------------------------
      integer function idKlaus(idPB)
      integer idPB
      idKlaus=0
      end function idKlaus


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for hqinitransfer
!-----------------------------------------------------------------------
      subroutine hqinitransfer
      end subroutine hqinitransfer

