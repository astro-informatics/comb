!------------------------------------------------------------------------------
! comb_tmplmap_mod -- COMB library real space template class
!
!! Contains definitions of template functions defined on the sky to initialise
!! comb obj objects.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   Septmeber 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module comb_tmplmap_mod

  use s2_types_mod
  use comb_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    comb_tmplmap_cos_thetaon2, &
    comb_tmplmap_point, &
    comb_tmplmap_butterfly, &
    comb_tmplmap_gaussian, &
    comb_tmplmap_mexhat, &
    comb_tmplmap_morlet


  !----------------------------------------------------------------------------

  contains

    
    !--------------------------------------------------------------------------
    ! comb_tmplmap_cos_thetaon2
    !
    !! Template function defined on the sphere.  
    !!   f(theta,phi) = cos(theta/2.0e0)
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplmap_cos_thetaon2(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val
      
      val = cos(theta/2.0e0)

    end function comb_tmplmap_cos_thetaon2


    !--------------------------------------------------------------------------
    ! comb_tmplmap_point
    !
    !! Template function defined on the sphere.  
    !! Point source centered at the north pole.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_tmplmap_point(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: theta_tol = 0.6e-1  !pi/20.0e0

      if(present(param)) then
         if(size(param) /= 1) then
          call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, 'comb_tmplmap_point')
        end if
        theta_tol = param(1)
      end if

      if(abs(theta) < theta_tol) then
         val = 1.0e0
      else
         val = 0.0e0
      end if

    end function comb_tmplmap_point


    !--------------------------------------------------------------------------
    ! comb_tmplmap_butterfly
    !
    !! Template function defined on the sphere.  
    !! Butterfly is a Gaussian in y direction and first derrivative of 
    !! Gaussian in x direction.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplmap_butterfly(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: N, arg

      if(present(param)) then
        if(size(param) /= 2) then
          call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, 'comb_tmplmap_butterfly')
        end if
        sigma_x = param(1)
        sigma_y = param(2)
      end if
      
      N = 1e0 / (2.0e0 * pi * sigma_x**3 * sigma_y)
      arg = -2.0e0 * (tan(theta/2.0e0))**2 * &
        (sigma_y**2 * (cos(phi))**2 + sigma_x**2 * (sin(phi))**2) &
        / (sigma_x**2 * sigma_y**2)

      val = - 4.0e0 * N * tan(theta/2.0e0) * cos(phi) / (1.0e0 + cos(theta)) &
       * exp(arg)

    end function comb_tmplmap_butterfly


    !--------------------------------------------------------------------------
    ! comb_tmplmap_gaussian
    !
    !! Template function defined on the sphere.  
    !! 2d Gaussian.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplmap_gaussian(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: N, arg

      if(present(param)) then
        if(size(param) /= 2) then
          call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, 'comb_tmplmap_gaussian')
        end if
        sigma_x = param(1)
        sigma_y = param(2)
      end if
      
      N = 1 / (2.0e0 * pi * sigma_x * sigma_y)
      arg = -2.0e0 * (tan(theta/2.0e0))**2 * &
        (sigma_y**2 * (cos(phi))**2 + sigma_x**2 * (sin(phi))**2) &
        / (sigma_x**2 * sigma_y**2)

      val =  N / (1.0e0 + cos(theta)) * exp(arg)

    end function comb_tmplmap_gaussian


    !--------------------------------------------------------------------------
    ! comb_tmplmap_mexhat
    !
    !! Template function defined on the sphere.  
    !! Mexican hat is negative of Laplacian of 2D Gaussian.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplmap_mexhat(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: N, arg, sigma_xony_sqrd

      if(present(param)) then
        if(size(param) /= 2) then
          call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, 'comb_tmplmap_mexhat')
        end if
        sigma_x = param(1)
        sigma_y = param(2)
      end if
      
      sigma_xony_sqrd = (sigma_x / sigma_y)**2

      N = 1 / (pi * sigma_x**3 * sigma_y**3)

      arg = -2.0e0 * (tan(theta/2.0e0))**2 * &
        (sigma_y**2 * (cos(phi))**2 + sigma_x**2 * (sin(phi))**2) &
        / (sigma_x**2 * sigma_y**2)

      val =  N / (1.0e0 + cos(theta)) * exp(arg) &
      * ( sigma_x**2 + sigma_y**2 - 4.0e0 * (tan(theta/2.0e0))**2 &
          * (cos(phi)**2 / sigma_xony_sqrd  + sin(phi)**2 * sigma_xony_sqrd) )

    end function comb_tmplmap_mexhat


    !--------------------------------------------------------------------------
    ! comb_tmplmap_morlet
    !
    !! Template function defined on the sphere.  
    !! Real Morlet wavelet.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as k0 of wave vector, where
    !!     k=(k0,0).
    !!   - Error occurs if param array has length greater than 1.
    !!   - If no parameter array is given then k0=10
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplmap_morlet(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: k0 = 10.0e0

      if(present(param)) then
          ! If only one parameter then taken as k0, where k=(k0,0).
         if(size(param) == 1) then
            k0 = param(1)
         else
            call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, &
              'comb_tmplmap_morlet')
         end if
      end if
      
      val = 2.0e0 / (1.0e0 + cos(theta)) &
            * cos(sqrt(2.0e0) * k0 * tan(theta/2.0e0) * cos(phi)) &
            * exp(-2.0e0*(tan(theta/2.0e0)**2))

    end function comb_tmplmap_morlet



end module comb_tmplmap_mod
