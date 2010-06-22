!------------------------------------------------------------------------------
! comb_obj_mod -- COMB library obj class
!
!! Provides functionality to support and manipulate a compact object defined 
!! on the sky.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module comb_obj_mod

  use comb_error_mod
  use s2_types_mod
  use s2_sky_mod
  use s2_pl_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------
  
  public :: &
    comb_obj_init, &
    comb_obj_free, &
    comb_obj_conv, &
    comb_obj_compute_alm, &
    comb_obj_write_sky, &
    comb_obj_get_init, &
    comb_obj_get_sky, &
    comb_obj_get_amplitude, &
    comb_obj_get_dilation, &
    comb_obj_get_alpha, &
    comb_obj_get_beta, &
    comb_obj_get_gamma, &
    comb_obj_get_name, &
    comb_obj_get_beam_status



  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface comb_obj_init
     module procedure &
       comb_obj_init_template, &
       comb_obj_init_mother,   &
       comb_obj_init_copy
  end interface
  

  !---------------------------------------
  ! Global variables
  !---------------------------------------
  
  ! None.
  

  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: comb_obj
     private
     logical :: init = .false.
     type(s2_sky) :: sky
     real(s2_sp) :: amplitude = 0.0e0
     real(s2_sp) :: dilation = 1.0e0
     real(s2_sp) :: alpha=0.0e0, beta=0.0e0, gamma=0.0e0
     character(len=S2_STRING_LEN) :: name = 'Not specified'
     logical :: beam_status = .false.
  end type comb_obj
  

  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! comb_obj_init_template
    !
    !! Initialise an obj from a template function defined over the sky.
    !!
    !! Notes:
    !!   - Dilation may be norm preserving or not.  Presently this option is
    !!     hard-coded and set to .false. (i.e. norm preserving scale *not* 
    !!     appliled).  Could be extended to an input option, although for 
    !!     this application may never require norm preserving dilation.
    !!   - Dilation may also be two dimensional (i.e. anisotropic).  
    !!     Currently not extended to this case, but very easy to do if 
    !!     required.
    !!
    !! Variables:
    !!   - template_fun: Template function to evaluate over the sky.
    !!   - nside: Healpix sky resolution.
    !!   - amplitude: Amplitude of the obj (i.e. value to scale the template 
    !!     function by).
    !!   - [pix_scheme_in]: Pixelisation scheme to create sky in (if not
    !!     specified then DEFAULT_PIS_SCHEME is used).
    !!   - [dilation]: Dilation to scale template by.
    !!   - [alpha]: Alpha Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - [beta]: Beta Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - [gamma]: Gamma Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - [name]: String specifying name of object.
    !!   - obj: Initialised obj.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function comb_obj_init_template(template_fun, nside, amplitude, &
      pix_scheme_in, dilation, alpha, beta, gamma, name) result(obj)

      real(s2_sp), intent(in) :: amplitude
      integer, intent(in) :: nside
      integer, intent(in), optional :: pix_scheme_in
      real(s2_sp), intent(in), optional :: dilation, alpha, beta, gamma
      character(len=*), intent(in), optional :: name
      type(comb_obj) :: obj
      interface 
         function template_fun(theta, phi, param) result(val)
           use s2_types_mod
           real(s2_sp), intent(in) :: theta, phi
           real(s2_sp), intent(in), optional :: param(:)
           real(s2_sp) :: val
         end function template_fun
      end interface
      
      integer, parameter :: DEFAUL_PIX_SCHEME = S2_SKY_RING
      integer :: pix_scheme

      ! Check object not already initialised.
      if(obj%init) then
        call comb_error(COMB_ERROR_INIT, 'comb_obj_init_template')
        return
      end if

      if(present(pix_scheme_in)) then
         pix_scheme = pix_scheme_in
      else
         pix_scheme = DEFAUL_PIX_SCHEME
      end if

      ! Initialise sky.
      obj%sky = s2_sky_init(template_fun, nside, pix_scheme)

      ! Set object attributes.
      obj%amplitude = amplitude
      if(present(name)) obj%name = name
      if(present(dilation)) obj%dilation = dilation
      if(present(alpha)) obj%alpha = alpha
      if(present(beta)) obj%beta = beta
      if(present(gamma)) obj%gamma = gamma

      ! Scale so have correct amplitude.
      call s2_sky_scale(obj%sky, obj%amplitude )

      ! Perform dilation if required.
      if(present(dilation)) then
         ! Only 1D dilation for now.
         call s2_sky_dilate(obj%sky, dilation, dilation, .false.)  
      end if

      ! Perform rotation if required.
      if(present(alpha) .and. present(beta) .and. present(gamma)) then
         call s2_sky_rotate(obj%sky, alpha, beta, gamma)
      end if

      ! Set object as initialised.
      obj%init = .true.

    end function comb_obj_init_template


    !--------------------------------------------------------------------------
    ! comb_obj_init_mother
    !
    !! Initialise obj from a mother obj.  The initialised obj is a scaled, 
    !! dilated and rotated version of the mother.
    !!
    !! Notes:
    !!   - Nside and pix_scheme same as mother.
    !!
    !! Variables:
    !!   - mother: Mother obj to `copy' when constructing this obj.
    !!   - amplitude: Amplitude of the obj (i.e. value to scale the template 
    !!     function by).
    !!   - dilation: Dilation to scale template by.
    !!   - alpha: Alpha Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - beta: beta Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - gamma: Gamma Euler rotation angle applied to rotate the template 
    !!     funciton defined on the sphere.
    !!   - [name]: String specifying name of object.
    !!   - obj: Initialised obj.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function comb_obj_init_mother(mother, amplitude, dilation, &
      alpha, beta, gamma, name) result(obj)

      type(comb_obj), intent(in) :: mother
      real(s2_sp), intent(in) :: amplitude
      real(s2_sp), intent(in) :: dilation, alpha, beta, gamma
      character(len=*), intent(in), optional :: name
      type(comb_obj) :: obj

      real(s2_sp), parameter :: TOL = 1e-4, NO_DILATION = 1.0e0, &
        ZERO_ANGLE = 0.0e0

!real(s2_sp) :: power

      ! Check object not already initialised.
      if(obj%init) then
        call comb_error(COMB_ERROR_INIT, 'comb_obj_init_mother')
        return
      end if

      ! Check mother object initialised.
      if(.not. mother%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_init_mother', &
          comment_add='Mother not initiliased')
        return
      end if

      ! To be valid mother object require no dilation and rotation.
      if(       abs(mother%dilation - NO_DILATION) > TOL &
           .or. abs(mother%alpha - ZERO_ANGLE) > TOL &
           .or. abs(mother%beta  - ZERO_ANGLE) > TOL &
           .or. abs(mother%gamma - ZERO_ANGLE) > TOL      ) then
         call comb_error(COMB_ERROR_OBJ_MOTHER_INVALID, 'comb_obj_init_mother')
         return
      end if

      ! Copy sky from mother.
      obj%sky = s2_sky_init(mother%sky)

      ! Set object name.
      if(present(name)) then
         obj%name = name
      else
         obj%name = mother%name
      end if

      ! Scale so have correct amplitude.
      if(abs(amplitude - mother%amplitude) > TOL) then
         call s2_sky_scale(obj%sky, amplitude / mother%amplitude )
         obj%amplitude = amplitude
      else
         obj%amplitude = mother%amplitude
      end if

      ! Dilate.
      ! Only 1D dilation for now.
      call s2_sky_dilate(obj%sky, dilation, dilation, .false.)  
      obj%dilation = dilation

!power = s2_sky_power_map(obj%sky)
!write(*,*) 'map power: ', power

      ! Rotate.
      call s2_sky_rotate(obj%sky, alpha, beta, gamma)
      obj%alpha = alpha
      obj%beta = beta
      obj%gamma = gamma

      ! Set object as initialised.
      obj%init = .true.

    end function comb_obj_init_mother


    !--------------------------------------------------------------------------
    ! comb_obj_init_copy
    !
    !! Initialise an obj as a copy of another obj.
    !!
    !! Variables:
    !!   - orig: Original obj to copy.
    !!   - copy: New object initialised as a copy of orig.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function comb_obj_init_copy(orig) result(copy)
      
      type(comb_obj), intent(in) :: orig
      type(comb_obj) :: copy

      ! Check original object initialised.
      if(.not. orig%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_init_copy')
      end if 

      ! Check copy object not already initialised.
      if(copy%init) then
        call comb_error(COMB_ERROR_INIT, 'comb_obj_init_copy')
        return
      end if

      ! Copy object attributes.

      copy%amplitude = orig%amplitude
      copy%dilation = orig%dilation
      copy%alpha = orig%alpha
      copy%beta = orig%beta
      copy%gamma = orig%gamma
      copy%name = orig%name
      copy%beam_status = orig%beam_status

      copy%sky = s2_sky_init(orig%sky)

      copy%init = .true.

    end function comb_obj_init_copy


    !--------------------------------------------------------------------------
    ! comb_obj_free
    !
    !! Free all data associated with an initialised obj and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - obj: Obj to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_obj_free(obj)

      type(comb_obj), intent(inout) :: obj

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_free')
      end if 

      ! Free sky.
      if(s2_sky_get_init(obj%sky)) call s2_sky_free(obj%sky)

      ! Set other parameters to default values.
      obj%amplitude = 0.0e0
      obj%dilation = 1.0e0
      obj%alpha = 0.0e0
      obj%beta = 0.0e0
      obj%gamma = 0.0e0
      obj%name = 'empty'
      obj%beam_status = .false.
      obj%init = .false.
      
    end subroutine comb_obj_free


    !--------------------------------------------------------------------------
    ! comb_obj_conv
    !
    !! Apply a beam by convolving it with the obj sky.
    !!
    !! Variables:
    !!   - obj: Obj object containing sky to apply beam to.
    !!   - beam: Beam to apply.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine comb_obj_conv(obj, beam)

      type(comb_obj), intent(inout) :: obj
      type(s2_pl), intent(in) :: beam

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_conv')
      end if 

      ! Convolve beam with wnoise sky.
      call s2_sky_conv(obj%sky, beam)

      ! Set beam status.
      obj%beam_status = .true.

    end subroutine comb_obj_conv


    !--------------------------------------------------------------------------
    ! comb_obj_compute_alm
    !
    !! Compute the alms of the obj sky at the lmax and mmax specified.
    !!
    !! Variables:
    !!   - obj: Obj object to compute alms of sky.
    !!   - lmax: Healpix spherical harmonic lmax.
    !!   - mmax: Healpix spherical harmonic mmax.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine comb_obj_compute_alm(obj, lmax, mmax)

      type(comb_obj), intent(inout) :: obj
      integer, intent(in) :: lmax, mmax

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_compute_alm')
      end if 

      ! Compute obj sky alms.
      call s2_sky_compute_alm(obj%sky, lmax, mmax)

    end subroutine comb_obj_compute_alm


    !--------------------------------------------------------------------------
    ! comb_obj_write_sky
    !
    !! Write an obj sky map to a fits file.
    !!
    !! Variables:
    !!   - obj: Obj containing sky to write to fits file.
    !!   - filename: Output fits file name.
    !!   o[comment]: Optional comment appended to header of fits file if 
    !!     present.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_obj_write_sky(obj, filename, comment)

      type(comb_obj), intent(in) :: obj
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

     ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_write_sky')
      end if 

      call s2_sky_write_map_file(obj%sky, filename, comment)

    end subroutine comb_obj_write_sky


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! comb_obj_get_init
    !
    !! Get init variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_init(obj) result(init)

      type(comb_obj), intent(in) :: obj
      logical :: init

      init = obj%init

    end function comb_obj_get_init


    !--------------------------------------------------------------------------
    ! comb_obj_get_sky
    !
    !! Get sky variable from the passed comb obj.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the obj std sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - sky: Object sky variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_sky(obj) result(sky)

      type(comb_obj), intent(in) :: obj
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_sky')
      end if 

      sky = s2_sky_init(obj%sky)

    end function comb_obj_get_sky


    !--------------------------------------------------------------------------
    ! comb_obj_get_amplitude
    !
    !! Get amplitude variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - amplitude: Object amplitude variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_amplitude(obj) result(amplitude)

      type(comb_obj), intent(in) :: obj
      real(s2_sp) :: amplitude

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_amplitude')
      end if 

      amplitude = obj%amplitude

    end function comb_obj_get_amplitude


    !--------------------------------------------------------------------------
    ! comb_obj_get_dilation
    !
    !! Get dilation variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - dilation: Object dilation variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_dilation(obj) result(dilation)

      type(comb_obj), intent(in) :: obj
      real(s2_sp) :: dilation

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_dilation')
      end if 

      dilation = obj%dilation

    end function comb_obj_get_dilation


    !--------------------------------------------------------------------------
    ! comb_obj_get_alpha
    !
    !! Get alpha variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - alpha: Object alpha variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_alpha(obj) result(alpha)

      type(comb_obj), intent(in) :: obj
      real(s2_sp) :: alpha

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_alpha')
      end if 

      alpha = obj%alpha

    end function comb_obj_get_alpha


    !--------------------------------------------------------------------------
    ! comb_obj_get_beta
    !
    !! Get beta variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - beta: Object beta variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_beta(obj) result(beta)

      type(comb_obj), intent(in) :: obj
      real(s2_sp) :: beta

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_beta')
      end if 

      beta = obj%beta

    end function comb_obj_get_beta


    !--------------------------------------------------------------------------
    ! comb_obj_get_gamma
    !
    !! Get gamma variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - gamma: Object gamma variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_gamma(obj) result(gamma)

      type(comb_obj), intent(in) :: obj
      real(s2_sp) :: gamma

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_gamma')
      end if 

      gamma = obj%gamma

    end function comb_obj_get_gamma


    !--------------------------------------------------------------------------
    ! comb_obj_get_name
    !
    !! Get name variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - name: Object name variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_name(obj) result(name)

      type(comb_obj), intent(in) :: obj
      character(len=S2_STRING_LEN) :: name

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_name')
      end if 

      name = obj%name

    end function comb_obj_get_name


    !--------------------------------------------------------------------------
    ! comb_obj_get_beam_status
    !
    !! Get beam_status variable from the passed comb obj.
    !!
    !! Variables:
    !!   - obj: COMB obj object to get the variable of.
    !!   - beam_status: Object beam_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function comb_obj_get_beam_status(obj) result(beam_status)

      type(comb_obj), intent(in) :: obj
      logical :: beam_status

      ! Check object initialised.
      if(.not. obj%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_obj_get_beam_status')
      end if 

      beam_status = obj%beam_status

    end function comb_obj_get_beam_status


end module comb_obj_mod
