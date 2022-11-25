!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine xCTEQ
!
!> @param[in] x
!> @param[in] q2
!> @param[in] is
!> @return 0.
!-----------------------------------------------------------------------
      double precision function xCTEQ(x,q2,is) 
      real x, q2
      integer is
      xCTEQ=0.
      end function xCTEQ     


!-----------------------------------------------------------------------
!> @brief
!> dummy interface for subroutine xCTEQ2
!
!> @param[in] idu1
!> @param[in] xx
!> @param[in] q2
!> @param[in] is
!> @param[in] idu2
!> @param[in] idu3
!> @return 0.
!-----------------------------------------------------------------------
      double precision function xCTEQ2(idu1,xx,q2,is,idu2,idu3)
      real xx, q2
      integer is, idu1, idu2, idu3
      xCTEQ2=0.
      end function xCTEQ2
