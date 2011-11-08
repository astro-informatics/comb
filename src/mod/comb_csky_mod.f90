!------------------------------------------------------------------------------
! comb_csky_mod -- COMB library csky class
!
!! Provides funcitonality to support a COMB csky object.  The csky object 
!! defines a full sky map consisting of embedded compact objects (COMB obj 
!! objects) and optional primodial cmb and noise realisations.  Output maps may
!! be produced contained the sum of all these components or each component map
!! may also be written individually.  An output parameter file may also be 
!! written to specify the compact object positions and other parameters, 
!! in addition to the cmb, beam and noise properties.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 August 2004
!
! Revisions:
!   August 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module comb_csky_mod

  use s2_types_mod
  use s2_error_mod
  use s2_sky_mod
  use s2_pl_mod
  use s2_cmb_mod
  use s2_wnoise_mod
  use comb_error_mod
  use comb_obj_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    comb_csky_init, &
    comb_csky_free, &
    comb_csky_write_param, &
    comb_csky_write_sky_obj, &
    comb_csky_write_sky_full, &
    comb_csky_write_sky_cmb, &
    comb_csky_write_sky_wnoise, &
    comb_csky_set_cmb, &
    comb_csky_set_wnoise, &
    comb_csky_get_init, &
    comb_csky_get_nobj, &
    comb_csky_get_nside, &
    comb_csky_get_obj, &
    comb_csky_get_sky_obj, &
    comb_csky_get_sky_full, &
    comb_csky_get_cmb, &
    comb_csky_get_wnoise, &
    comb_csky_get_cmb_status, &
    comb_csky_get_wnoise_status


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface comb_csky_init
     module procedure &
       comb_csky_init_mother, &
       comb_csky_init_array
  end interface

  
  !---------------------------------------
  ! Global variables
  !---------------------------------------

  ! None.

  
  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: comb_csky
     private
     logical :: init = .false.
     integer :: nobj = 0
     integer :: nside = 0
     type(comb_obj), allocatable :: obj(:)
     type(s2_sky) :: sky_obj
     type(s2_sky) :: sky_full
     type(s2_cmb) :: cmb
     type(s2_wnoise) :: wnoise
     logical :: cmb_status = .false.
     logical :: wnoise_status = .false.
  end type comb_csky


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! comb_csky_init_mother
    !
    !! Initialise a csky from a mother obj and arrays of parameters specifying 
    !! properties of each object to add to the obj array.  Cmb and wnoise 
    !! objects may also be optionally added.
    !!
    !! Variables:
    !!   - obj: Mother obj to generate other objs from.
    !!   - nobj: Number objects to be created (must be same as size of all 
    !!     parameter arrays).
    !!   - amplitude: Array of amplitudes for objs constructed from mother.
    !!   - dilation: Array of dilations for objs constructed from mother.
    !!   - alpha: Array of Euler alpha angles for objs constructed from mother.
    !!   - beta: Array of Euler beta angles for objs constructed from mother.
    !!   - gamma: Array of Euler gamma angles for objs constructed from mother.
    !!   - [cmb]: Cmb to add to csky if present.
    !!   - [wnoise]: Wnoise to add to csky if present.
    !!   - [beam]: Beam to apply to individual objects.  If beam is specified 
    !!     lmax and mmax must also be specified.
    !!   - [lmax]: Healpix lmax for computing alms of objects when applying 
    !!     beam.
    !!   - [mmax]: Healpix mmax for computing alms of objects when applying 
    !!     beam.
    !!   - csky: Initialised csky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function comb_csky_init_mother(obj, nobj, amplitude, dilation, &
         alpha, beta, gamma, cmb, wnoise, beam, lmax, mmax) result(csky)

      type(comb_obj), intent(in) :: obj
      integer, intent(in) :: nobj
      real(s2_sp), intent(in) :: amplitude(:), dilation(:)
      real(s2_sp), intent(in) :: alpha(:), beta(:), gamma(:)
      type(s2_cmb), intent(in), optional :: cmb
      type(s2_wnoise), intent(in), optional :: wnoise
      type(s2_pl), intent(in), optional :: beam
      integer, intent(in), optional :: lmax, mmax
      type(comb_csky) :: csky

      integer :: iobj, fail
      type(s2_sky) :: temp_sky

      ! Check object not already initialised.
      if(csky%init) then
        call comb_error(COMB_ERROR_INIT, 'comb_csky_init_mother')
        return
      end if

      ! Check arrays have same length as nobj.
      if( size(amplitude) /= nobj .or. &
          size(dilation)  /= nobj .or. &
          size(alpha)     /= nobj .or. &
          size(beta)      /= nobj .or. &
          size(gamma)     /= nobj ) then
         call comb_error(COMB_ERROR_INIT_FAIL, 'comb_sky_init_mother', &
           comment_add='Inconsistent size for parameter arrays')
      end if

      ! Initialise attributes.
      csky%nobj = nobj 
      temp_sky = comb_obj_get_sky(obj)
      csky%nside = s2_sky_get_nside(temp_sky)
      call s2_sky_free(temp_sky)

      ! Initialise compact objects.
      allocate(csky%obj(1:csky%nobj), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_csky_init_mother')
      end if
      do iobj = 1,csky%nobj
         csky%obj(iobj) = comb_obj_init(obj, amplitude(iobj), dilation(iobj), &
              alpha(iobj), beta(iobj), gamma(iobj))
         if(present(beam)) then
            if(present(lmax) .and. present(mmax)) then
               call comb_obj_compute_alm(csky%obj(iobj), lmax, mmax)
               call comb_obj_conv(csky%obj(iobj), beam)
            else
               call comb_error(COMB_ERROR_CSKY_LMAX_NOT_DEF, &
                 'comb_sky_init_mother')
            end if
         end if
      end do

      ! Set object as initialised.
      csky%init = .true.

      ! Make obj sky.
      call comb_csky_compute_sky_obj(csky)

      ! Calculate full sky.
      call comb_csky_compute_sky_full(csky)

      ! Set cmb if present
      if(present(cmb)) then
         call comb_csky_set_cmb(csky, cmb)
      end if

      ! Set wnoise if present.
      if(present(wnoise)) then
         call comb_csky_set_wnoise(csky, wnoise) 
      end if

    end function comb_csky_init_mother


    !--------------------------------------------------------------------------
    ! comb_csky_init_array
    !
    !! Initialise csky from an array of already defined objs.  Cmb and wnoise 
    !! objects may also be optionally added.
    !!
    !! Variables:
    !!   - obj: Array of obj objects to add to csky.
    !!   - [cmb]: Cmb to add to csky if present.
    !!   - [wnoise]: Wnoise to add to csky if present.
    !!   - csky: Initialised csky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function comb_csky_init_array(obj, cmb, wnoise, beam, &
         lmax, mmax) result(csky)

      type(comb_obj), intent(in) :: obj(:)
      type(s2_cmb), intent(in), optional :: cmb
      type(s2_wnoise), intent(in), optional :: wnoise
      type(s2_pl), intent(in), optional :: beam
      integer, intent(in), optional :: lmax, mmax
      type(comb_csky) :: csky

      integer :: iobj, fail
      integer :: nside, pix_scheme
      integer :: lmax_test, mmax_test
      logical :: harmonic_tmpl
      type(s2_sky) :: temp_sky

      ! Check object not already initialised.
      if(csky%init) then
        call comb_error(COMB_ERROR_INIT, 'comb_csky_init_array')
        return
      end if

      ! Set attributes.

      csky%nobj = size(obj)
      if(csky%nobj <= 0 ) then
         call comb_error(COMB_ERROR_INIT_FAIL, 'comb_csky_init_array', &
              comment_add='nobj <= 0')
      end if
      
      temp_sky = comb_obj_get_sky(obj(1))
      csky%nside = s2_sky_get_nside(temp_sky)
      call s2_sky_free(temp_sky)

      ! Initialise compact objects.

      allocate(csky%obj(1:csky%nobj), stat=fail)
      if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_csky_init_array')
      end if

      ! Get properties of first object sky to later check all objects skies
      ! are consistent.
      temp_sky = comb_obj_get_sky(obj(1))
      nside = s2_sky_get_nside(temp_sky)
      if(present(lmax)) then
         lmax_test = lmax
      else
         lmax_test = s2_sky_get_lmax(temp_sky)
      end if
      if(present(mmax)) then
         mmax_test = mmax
      else
         mmax_test = s2_sky_get_mmax(temp_sky)
      end if
      pix_scheme = s2_sky_get_pix_scheme(temp_sky)
      harmonic_tmpl = comb_obj_get_harmonic_tmpl(obj(1))
      call s2_sky_free(temp_sky)

      do iobj = 1,csky%nobj
         ! Check each object has same properties.
         temp_sky = comb_obj_get_sky(obj(iobj))
         if(s2_sky_get_pix_scheme(temp_sky) /= pix_scheme &
            .or. s2_sky_get_nside(temp_sky) /= nside & 
            .or. s2_sky_get_lmax(temp_sky) /= lmax_test & 
            .or. s2_sky_get_mmax(temp_sky) /= mmax_test &
            .or. comb_obj_get_harmonic_tmpl(obj(iobj)) .neqv. harmonic_tmpl ) then 
           call comb_error(COMB_ERROR_INIT_FAIL, 'comb_csky_init_array', &
             comment_add=&
             'Compact objects in array have inconsistent properties')
         end if
         call s2_sky_free(temp_sky)
         ! Copy compact object.
         csky%obj(iobj) = comb_obj_init(obj(iobj))
         
         ! Apply beam if present.
         if(present(beam)) then
            call comb_obj_compute_alm(csky%obj(iobj), lmax_test, mmax_test)
            call comb_obj_conv(csky%obj(iobj), beam)
         end if

      end do

      ! Set object as initialised.
      csky%init = .true.

      ! Make obj sky.
      call comb_csky_compute_sky_obj(csky)

      ! Calculate full sky.
      call comb_csky_compute_sky_full(csky)

      ! Set cmb if present.
      if(present(cmb)) then
         call comb_csky_set_cmb(csky, cmb) 
      end if
      
      ! Set wnoise if present.
      if(present(wnoise)) then
         call comb_csky_set_wnoise(csky, wnoise) 
      end if
      
    end function comb_csky_init_array


    !--------------------------------------------------------------------------
    ! comb_csky_free
    !
    !! Free all data associated with an initialised csky and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - csky: Csky to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_free(csky)
    
      type(comb_csky), intent(inout) :: csky
      
      integer :: iobj

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_free')
      end if 
      
      ! Free compact objects.
      do iobj = 1,csky%nobj
         call comb_obj_free(csky%obj(iobj))
      end do
      deallocate(csky%obj)

      ! Free skies.
      if(s2_sky_get_init(csky%sky_obj)) call s2_sky_free(csky%sky_obj)
      if(s2_sky_get_init(csky%sky_full)) call s2_sky_free(csky%sky_full)

      ! Free cmb and wnoise if they exist.
      if(csky%cmb_status) call s2_cmb_free(csky%cmb)
      if(csky%wnoise_status) call s2_wnoise_free(csky%wnoise)

      ! Reset other object attributes.
      csky%nobj = 0
      csky%nside = 0
      csky%cmb_status = .false.
      csky%wnoise_status = .true.
      csky%init = .false.

    end subroutine comb_csky_free


    !--------------------------------------------------------------------------
    ! comb_csky_compute_sky_obj
    !
    !! Compute the full sky obj representation consisting of the sum of all
    !! individual objs contained in the csky obj array.
    !!
    !! Variables:
    !!   - csky: Csky to compute full sky obj representation of.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_compute_sky_obj(csky)

      type(comb_csky), intent(inout) :: csky

      integer :: iobj
      type(s2_sky) :: obj_sky_temp, csky_sky_temp
      logical :: harmonic_tmpl

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_free')
      end if 

      ! Initialise object sky with first object.
      obj_sky_temp = comb_obj_get_sky(csky%obj(1))
      csky%sky_obj = s2_sky_init(obj_sky_temp)
      call s2_sky_free(obj_sky_temp)

      ! Then add all other compact objects to sky.
      do iobj = 2,csky%nobj
         obj_sky_temp = comb_obj_get_sky(csky%obj(iobj))
         harmonic_tmpl = comb_obj_get_harmonic_tmpl(csky%obj(iobj))
         if (harmonic_tmpl) then
            csky_sky_temp = s2_sky_add_alm(csky%sky_obj, obj_sky_temp)
         else
            csky_sky_temp = s2_sky_add(csky%sky_obj, obj_sky_temp)
         end if
         call s2_sky_free(csky%sky_obj)
         csky%sky_obj = s2_sky_init(csky_sky_temp)
         call s2_sky_free(csky_sky_temp)
         call s2_sky_free(obj_sky_temp)
      end do

    end subroutine comb_csky_compute_sky_obj


    !--------------------------------------------------------------------------
    ! comb_csky_compute_sky_full
    !
    !! Compute full sky map consisting of the full sky obj representation and 
    !! the cmb and noise realisation (if present).
    !!
    !! Variables:
    !!   - csky: Csky to compute full sky map (sky_full) consisting of all
    !!     present components.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_compute_sky_full(csky)

      type(comb_csky), intent(inout) :: csky

      type(s2_sky) :: sky_temp, csky_sky_temp

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_free')
      end if 

      ! If sky_full already initialised then free.
      if(s2_sky_get_init(csky%sky_full)) call s2_sky_free(csky%sky_full)
      
      ! Initialise as copy of compact object sky, then add other skies
      ! if present.
      csky%sky_full = s2_sky_init(csky%sky_obj)

      ! If cmb present then add cmb sky.
      if(csky%cmb_status) then
         sky_temp = s2_cmb_get_sky(csky%cmb)
         csky_sky_temp = s2_sky_add(csky%sky_full, sky_temp)
         call s2_sky_free(csky%sky_full)
         csky%sky_full = s2_sky_init(csky_sky_temp)
         call s2_sky_free(csky_sky_temp)
         call s2_sky_free(sky_temp)
      end if

      ! If wnoise present then add wnoise sky.
      if(csky%wnoise_status) then
         sky_temp = s2_wnoise_get_sky(csky%wnoise)
         csky_sky_temp = s2_sky_add(csky%sky_full, sky_temp)
         call s2_sky_free(csky%sky_full)
         csky%sky_full = s2_sky_init(csky_sky_temp)
         call s2_sky_free(csky_sky_temp)
         call s2_sky_free(sky_temp)
      end if

    end subroutine comb_csky_compute_sky_full


    !--------------------------------------------------------------------------
    ! comb_csky_write_param
    !
    !! Write text parameter file describing csky attributes.  All compact 
    !! object parameters are written (eg. alplitude, dilation and Euler 
    !! angles), plus limited cmb and noise parameters).
    !!
    !! Variables:
    !!   - csky: Csky contained parameters to be written.
    !!   - filename: Name of output parameter file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_write_param(csky, filename)

      type(comb_csky), intent(in) :: csky
      character(len=*), intent(in) :: filename

      character(len=8) :: date
      character(len=10) :: time
      integer :: fileid = 10, iobj

      ! Open file.
      open(unit=fileid, file=filename, status='replace', action='write')

      !---------------------------------------
      ! Write file header
      !---------------------------------------

      write(fileid, '(a)') '# COMB_CSKY Output Parameter Summary'
      call date_and_time(date, time)
      write(fileid, '(a,a,a,a,a,a,a,a,a,a,a,a)') '# Produced: ', &
           time(1:2), ':', time(3:4), ':', time(5:6), ' on ', &
           date(7:8), '/', date(5:6), '/', date(1:4)
      write(fileid, '(a)') '# ----------------------------------'
      write(fileid, *)

      !---------------------------------------
      ! Write healpix parameters
      !---------------------------------------

      write(fileid, '(a)') '# Healpix parameters'
      write(fileid, '(a)') '# ----------------------------------'
      write(fileid, '(a,i13)') 'nside = ', csky%nside
      if(csky%cmb_status) then
         write(fileid, '(a,i14)') 'lmax = ', &
           s2_cmb_get_lmax(csky%cmb)
      end if
      write(fileid, *)

      !---------------------------------------
      ! Write obj parameters
      !---------------------------------------

      write(fileid, '(a)') '# Compact object parameters'
      write(fileid, '(a)') '# ----------------------------------'

      do iobj = 1, csky%nobj
 
         write(fileid, '(a,i3,a)') 'object', iobj, ':'
         write(fileid, '(a,a13)') ' name = ', &
           trim(comb_obj_get_name(csky%obj(iobj)))
         write(fileid, '(a,f8.4)') ' amplitude = ', &
           comb_obj_get_amplitude(csky%obj(iobj))
         write(fileid, '(a,f9.4)') ' dilation = ', &
           comb_obj_get_dilation(csky%obj(iobj))
         write(fileid, '(a,f12.4)') ' alpha = ', &
           comb_obj_get_alpha(csky%obj(iobj))
         write(fileid, '(a,f13.4)') ' beta = ', &
           comb_obj_get_beta(csky%obj(iobj))
         write(fileid, '(a,f12.4)') ' gamma = ', &
           comb_obj_get_gamma(csky%obj(iobj))

      end do
      write(fileid, *)

      !---------------------------------------
      ! Write cmb parameters
      !---------------------------------------

      write(fileid, '(a)') '# CMB parameters'
      write(fileid, '(a)') '# ----------------------------------'
      write(fileid, '(a,l8)') 'cmb_status = ', csky%cmb_status
      write(fileid, *)

      !---------------------------------------
      ! Write beam parameters
      !---------------------------------------

      write(fileid, '(a)') '# Beam parameters'
      write(fileid, '(a)') '# ----------------------------------'
      if(csky%cmb_status) then
         write(fileid, '(a,l7)') 'beam_status = ', &
              s2_cmb_get_beam_applied(csky%cmb)
      else
         write(fileid, '(a,l7)') 'beam_status = ', .false.
      end if
      write(fileid, *)

      !---------------------------------------
      ! Write noise parameters
      !---------------------------------------

      write(fileid, '(a)') '# White noise parameters'
      write(fileid, '(a)') '# ----------------------------------'
      write(fileid, '(a,l6)') 'noise_status = ', csky%wnoise_status
      if(csky%wnoise_status) then
         if(s2_wnoise_get_type(csky%wnoise) == S2_WNOISE_TYPE_STD_SKY) then
            write(fileid, '(a,l4)') 'noise_full_sky = ', .true.
            write(fileid, '(a,f6.4)') 'noise_sigma0 = ', &
              s2_wnoise_get_sigma0(csky%wnoise)
         else
            write(fileid, '(a,l4)') 'noise_full_sky = ', .false.
            write(fileid, '(a,f9.4)') 'noise_std = ', &
              s2_wnoise_get_std_const(csky%wnoise)
         end if      
      end if

      write(fileid, *)

      ! Close file.
      close(unit=fileid)

    end subroutine comb_csky_write_param


    !--------------------------------------------------------------------------
    ! comb_csky_write_sky_obj
    !
    !! Write csky compact object full sky representation to fits file.
    !!
    !! Variables:
    !!   - csky: Csky containing sky_obj to be written to file.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !!   - [file_type_in]: Optional specifying type of file to write.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !   June 2010 - Ability to write different file types added by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_write_sky_obj(csky, filename, comment, file_type_in)

      type(comb_csky), intent(inout) :: csky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment
      integer, intent(in), optional :: file_type_in

      integer :: file_type = S2_SKY_FILE_TYPE_MAP

      if(present(file_type_in)) file_type = file_type_in
      
      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_obj')
      end if 

      ! Ensure sky computed in required space (map or alms).
      select case (file_type)
        case(S2_SKY_FILE_TYPE_MAP)
           call s2_sky_compute_map(csky%sky_obj)
        case(S2_SKY_FILE_TYPE_ALM)
           call s2_sky_compute_alm(csky%sky_obj)
      end select

      ! Write file.
      call s2_sky_write_file(csky%sky_obj, filename, file_type, comment)

    end subroutine comb_csky_write_sky_obj


    !--------------------------------------------------------------------------
    ! comb_csky_write_sky_full
    !
    !! Write csky full sky representation consisting of all present components
    !! to fits file.
    !!
    !! Variables:
    !!   - csky: Csky containing sky_full to be written to output file.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !!   - [file_type_in]: Optional specifying type of file to write.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !   June 2010 - Ability to write different file types added by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_write_sky_full(csky, filename, comment, file_type_in)

      type(comb_csky), intent(inout) :: csky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment
      integer, intent(in), optional :: file_type_in

      integer :: file_type = S2_SKY_FILE_TYPE_MAP

      if(present(file_type_in)) file_type = file_type_in

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_full')
      end if 

      ! Ensure sky computed in required space (map or alms).
      select case (file_type)
        case(S2_SKY_FILE_TYPE_MAP)
           call s2_sky_compute_map(csky%sky_full)
        case(S2_SKY_FILE_TYPE_ALM)
           call s2_sky_compute_alm(csky%sky_full)
      end select

      ! Write file.
      call s2_sky_write_file(csky%sky_full, filename, file_type, comment)

    end subroutine comb_csky_write_sky_full


    !--------------------------------------------------------------------------
    ! comb_csky_write_sky_cmb
    !
    !! Write cmb sky map (if present) to fits file.
    !!
    !! Variables:
    !!   - csky: Csky containing the cmb, that in turn contains the sky 
    !!     realisation to be written.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_write_sky_cmb(csky, filename, comment)

      type(comb_csky), intent(in) :: csky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_cmb')
      end if 

      if(csky%cmb_status) then
         call s2_cmb_write_sky(csky%cmb, filename, comment)
      else
         call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_cmb', &
           comment_add='Warning: cmb not present so not written', &
           halt_in=.false.)
      end if

    end subroutine comb_csky_write_sky_cmb


    !--------------------------------------------------------------------------
    ! comb_csky_write_sky_wnoise
    !
    !! Write wnoise sky map (if present) to fits file.
    !!
    !! Variables: 
    !!   - csky: Csky containing the wnoise, that in turn contains the sky
    !!     realisation to be written.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits file
    !!     header.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_write_sky_wnoise(csky, filename, comment)

      type(comb_csky), intent(in) :: csky
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_wnoise')
      end if 

      if(csky%wnoise_status) then
         call s2_wnoise_write_sky_file(csky%wnoise, filename, comment)
      else
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_write_sky_wnoise', &
          comment_add='Warning: wnoise not present so not written', &
          halt_in=.false.)
      end if

    end subroutine comb_csky_write_sky_wnoise


    !--------------------------------------------------------------------------
    ! Set routines
    !--------------------------------------------------------------------------
   
    !--------------------------------------------------------------------------
    ! comb_csky_set_cmb
    !
    !! Set the csky cmb object to that passed as an argument (actually 
    !! create a copy of this cmb).  The full sky csky map is then recomputed 
    !! with the new cmb included.
    !!
    !! Variables:
    !!   - csky: Csky to set cmb of.
    !!   - cmb: New csky cmb to be set.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_set_cmb(csky, cmb)
      
      type(comb_csky), intent(inout) :: csky
      type(s2_cmb), intent(in) :: cmb

      ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_set_cmb')
      end if 

      ! If already initialised then free.
      if(s2_cmb_get_init(csky%cmb)) call s2_cmb_free(csky%cmb)

      ! Replace with new passed cmb.
      csky%cmb = s2_cmb_init(cmb)

      csky%cmb_status = .true.

      ! Re-compute full sky.
      call comb_csky_compute_sky_full(csky)

    end subroutine comb_csky_set_cmb


    !--------------------------------------------------------------------------
    ! comb_csky_set_wnoise
    !
    !! Set the csky wnoise object to that passed as an argument (actually 
    !! create a copy of this wnoise).  The full sky csky map is then
    !! recomputed with the new wnoise included.
    !!
    !! Variables:
    !!   - csky: Csky to set wnoise of.
    !!   - wnoise: New csky wnoise to be set.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine comb_csky_set_wnoise(csky, wnoise)

      type(comb_csky), intent(inout) :: csky
      type(s2_wnoise) :: wnoise
      
     ! Check object initialised.
      if(.not. csky%init) then
        call comb_error(COMB_ERROR_NOT_INIT, 'comb_csky_set_wnoise')
      end if 

      ! If already initialised then free.
      if(s2_wnoise_get_init(csky%wnoise)) call s2_wnoise_free(csky%wnoise)

      ! Replace with new passed wnoise.
      csky%wnoise = s2_wnoise_init(wnoise)

      csky%wnoise_status = .true.

      ! Re-compute full sky.
      call comb_csky_compute_sky_full(csky)

    end subroutine comb_csky_set_wnoise


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------
   
    !--------------------------------------------------------------------------
    ! comb_csky_get_init
    !
    !! Get init variable from the passed csky object.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - init: Object init variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_init(csky) result(init)
      
      type(comb_csky), intent(in) :: csky
      logical :: init

      init = csky%init

    end function comb_csky_get_init


    !--------------------------------------------------------------------------
    ! comb_csky_get_nobj
    !
    !! Get nobj variable from the passed csky object.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - nobj: Object nobj variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_nobj(csky) result(nobj)
      
      type(comb_csky), intent(in) :: csky
      integer :: nobj

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_nobj')
      end if

      nobj = csky%nobj

    end function comb_csky_get_nobj


    !--------------------------------------------------------------------------
    ! comb_csky_get_nside
    !
    !! Get nside variable from the passed csky object.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - nside: Object nside variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_nside(csky) result(nside)
      
      type(comb_csky), intent(in) :: csky
      integer :: nside

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_nside')
      end if

      nside = csky%nside

    end function comb_csky_get_nside


    !--------------------------------------------------------------------------
    ! comb_csky_get_obj
    !
    !! Get compact object variable at specified index from the passed 
    !! csky object.
    !!    
    !! Notes:
    !!   - Initialises a new obj as a copy of the csky obj.
    !!   - The returned std obj is subsequently independed of the obj stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - iobj: Index of compact object to get.
    !!   - obj: Object compact object (obj) object returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_obj(csky, iobj) result(obj)
      
      type(comb_csky), intent(in) :: csky
      integer, intent(in) :: iobj
      type(comb_obj) :: obj

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_obj')
      end if

      obj = comb_obj_init(csky%obj(iobj))

    end function comb_csky_get_obj


    !--------------------------------------------------------------------------
    ! comb_csky_get_sky_obj
    !
    !! Get sky_obj variable from the passed csky object.
    !! 
    !! Notes:
    !!   - Initialises a new sky as a copy of the obj sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - sky_obj: Object sky_obj variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_sky_obj(csky) result(sky_obj)
      
      type(comb_csky), intent(in) :: csky
      type(s2_sky) :: sky_obj

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_sky_obj')
      end if

      sky_obj = s2_sky_init(csky%sky_obj)

    end function comb_csky_get_sky_obj


    !--------------------------------------------------------------------------
    ! comb_csky_get_sky_full
    !
    !! Get sky_full variable from the passed csky object.
    !!   
    !! Notes:
    !!   - Initialises a new sky as a copy of the csky full sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - sky_full: Object sky_full variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_sky_full(csky) result(sky_full)
      
      type(comb_csky), intent(in) :: csky
      type(s2_sky) :: sky_full

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_sky_full')
      end if

      sky_full = s2_sky_init(csky%sky_full)

    end function comb_csky_get_sky_full


    !--------------------------------------------------------------------------
    ! comb_csky_get_cmb
    !
    !! Get cmb variable from the passed csky object.
    !!
    !! Notes:
    !!   - Initialises a new cmb as a copy of the csky cmb.
    !!   - The returned cmb is subsequently independed of the cmb stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - cmb: Object cmb variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_cmb(csky) result(cmb)
      
      type(comb_csky), intent(in) :: csky
      type(s2_cmb) :: cmb

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_cmb')
      end if

      cmb = s2_cmb_init(csky%cmb)

    end function comb_csky_get_cmb


    !--------------------------------------------------------------------------
    ! comb_csky_get_wnoise
    !
    !! Get wnoise variable from the passed csky object.
    !!
    !! Notes:
    !!   - Initialises a new wnoise as a copy of the csky wnoise.
    !!   - The returned wnoise is subsequently independed of the wnoise stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - wnoise: Object wnoise variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_wnoise(csky) result(wnoise)
      
      type(comb_csky), intent(in) :: csky
      type(s2_wnoise) :: wnoise

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_wnoise')
      end if

      wnoise = s2_wnoise_init(csky%wnoise)

    end function comb_csky_get_wnoise


    !--------------------------------------------------------------------------
    ! comb_csky_get_cmb_status
    !
    !! Get cmb_status variable from the passed csky object.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - cmb_status: Object cmb_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_cmb_status(csky) result(cmb_status)
      
      type(comb_csky), intent(in) :: csky
      logical :: cmb_status

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_cmb_status')
      end if

      cmb_status = csky%cmb_status

    end function comb_csky_get_cmb_status


    !--------------------------------------------------------------------------
    ! comb_csky_get_wnoise_status
    !
    !! Get wnoise_status variable from the passed csky object.
    !!
    !! Variables:
    !!   - csky: COMB csky object to get the variable of.
    !!   - wnoise_status: Object wnoise_status variable returned.
    !
    !! @author J. D. McEwen
    !! @version 0.1 August 2004
    !
    ! Revisions:
    !   August 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function comb_csky_get_wnoise_status(csky) result(wnoise_status)
      
      type(comb_csky), intent(in) :: csky
      logical :: wnoise_status

      ! Check object initialised.
      if(.not. csky%init) then
        call s2_error(COMB_ERROR_NOT_INIT, 'comb_csky_get_wnoise_status')
      end if

      wnoise_status = csky%wnoise_status

    end function comb_csky_get_wnoise_status


end module comb_csky_mod
