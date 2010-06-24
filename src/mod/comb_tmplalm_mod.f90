!------------------------------------------------------------------------------
! comb_tmpl_mod -- COMB library template class
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

module comb_tmplalm_mod

  use s2_types_mod
  use comb_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    comb_tmplalm_gaussian


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! comb_tmplalm_gaussian
    !
    !! Template function defined on the sphere in harmonic space.  
    !! 2d Gaussian.
    !!
    !! Variables:
    !!   - el: Spherical harmonic order l.
    !!   - m: Sperhical harmonic order m.
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version Under svn version control
    !
    ! Revisions:
    !   June 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function comb_tmplalm_gaussian(el, m, param) result(val)

      integer, intent(in) :: el, m
      real(s2_sp), intent(in), optional :: param(:)
      complex(s2_spc) :: val

      real(s2_sp) :: sigma = 0.05e0
      real(s2_sp) :: N, arg

      if(present(param)) then
        if(size(param) /= 1) then
          call comb_error(COMB_ERROR_TMPL_PARAM_INVALID, 'comb_tmplalm_gaussian')
        end if
        sigma = param(1)
      end if

      if (m == 0) then
         val =  exp(-el*el*sigma*sigma/2e0)
      else
         val = 0
      end if

    end function comb_tmplalm_gaussian



end module comb_tmplalm_mod
