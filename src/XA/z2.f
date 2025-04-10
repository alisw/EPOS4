!--------------------------------------------------------------------------
!
!    This piece of code is meant to be changed by the user
!
!--------------------------------------------------------------------------

!-----------------------------------------------------------------------
! The following example subroutine gives access to the hydro table 
! It is called after hydro, before hadronization, and before hadronic cascade)
!     
!    attention:
!
!        don't move "call pitabget(...)" outside this subroutine !!
!            The corresponding C++ object is destroyed after  the 
!                call to useHydroTable
!-----------------------------------------------------------------------
      subroutine useHydroTable
#include "aaa.h"
#include "ho.h"
      double precision temc
      common/ctemc/temc(netahxx,ntauhxx,nxhxx,nyhxx)
      common/cpssc/pssc(netahxx,ntauhxx,nxhxx,nyhxx)
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! dimensions: 
      !
      !     ntauhy,nxhy,nyhy,nzhy  (z refers to eta)
      !
      ! ranges: 
      !
      !     tauminhy,taumaxhy,xminhy,xmaxhy,yminhy,ymaxhy,zminhy,zmaxhy
      ! 
      ! for given ntau,nx,ny,nz the following table elements may be used:
      !
      !   epsc(nz,ntau,nx,ny)      (energy density)
      !   temc(nz,ntau,nx,ny)      (temperature)
      !   pssc(nz,ntau,nx,ny)      (pressure)
      !   velc(1,nz,ntau,nx,ny)    (velocity x comp)
      !   velc(2,nz,ntau,nx,ny)    (velocity y comp)
      !   velc(3,nz,ntau,nx,ny)    (velocity z comp)
      !   barc(1,nz,ntau,nx,ny)    (baryon density)
      !   barc(2,nz,ntau,nx,ny)    (charge density)
      !   barc(3,nz,ntau,nx,ny)    (strange density)
      !
      ! pi_mu_nu (shear stress tensor)  may be obtained via
      !
      !   call pitabget(mu,nu,nz,ntau,nx,ny, pi_mu_nu )
      !
      !  mu and nu refer to: 1=tau, 2=x, 3=y, 4=z(=eta) 
      !
      ! pi_mu_nu are components in Milne coordinates,
      ! with each eta component muliplied by tau,
      ! so for example pi_mu_nu is
      !           pi^{tau tau}  for mu=1,nu=1
      !           pi^{tau  x }  for mu=1,nu=2
      !    tau  * pi^{tau eta}  for mu=1,nu=4
      !  tau**2 * pi^{eta eta}  for mu=4,nu=4
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      
      !    .... put your analysis here ....
      
      
      end

