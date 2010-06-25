!------------------------------------------------------------------------------
! comb_csim -- COMB program
!
!! Generates simulations of compact objects embedded in primordial cmb 
!! (convolved with beam if desired) and noise.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

program comb_csim

  use s2_types_mod, only: S2_STRING_LEN, s2_sp, pi
  use s2_error_mod
  use s2_distn_mod
  use s2_cmb_mod
  use s2_pl_mod
  use s2_vect_mod, only: s2_vect_arcmin_to_rad
  use s2_wnoise_mod
  use s2_sky_mod
  use comb_error_mod
  use comb_tmplmap_mod
  use comb_tmplalm_mod
  use comb_obj_mod
  use comb_csky_mod
  
  use extension, only: getArgument, nArguments
  use paramfile_io, only: paramfile_handle, parse_init, parse_int, &
    parse_real, parse_lgt, parse_string, concatnl
  use pix_tools, only: nside2npix

  implicit none

  !----------------------------------------------------------------------------

  character(len=*), parameter :: CODE_NAME = 'COMB_CSIM'
  character(len=*), parameter :: CODE_DESCRIPTION = &
    'Compact object embedded CMB simulation'
  character(len=*), parameter :: CODE_VERSION = '1.0'
  character(len=*), parameter :: CODE_AUTHOR = 'Jason McEwen'

  character(len=S2_STRING_LEN) :: filename_param
  character(len=S2_STRING_LEN) :: description
  character(len=S2_STRING_LEN) :: line, line2
  type(paramfile_handle) :: handle

  integer, parameter :: NUM_COMMENT_LINES_CLT_FILE = 29   ! Value for WMAP1
!  integer, parameter :: NUM_COMMENT_LINES_CLT_FILE = 45   ! Value for WMAP3
  integer, parameter :: NUM_COMMENT_LINES_BEAM_FILE = 8 

  real(s2_sp) :: TOL = 1e-5

  !----------------------------------------------------------------------------

  ! Define default and allowable parameters.

  ! If integer or real then OPT (option) array specifies lower and upper
  ! bounds allowed.
  ! If string then OPT (option) array specified allowable strings.

  !---------------------------------------

  ! Healpix parameters.

  integer :: nside, lmax, mmax
  integer, parameter :: lmin = 2  ! Monopole and dilpole excluded.

  integer, parameter ::     &
    NSIDE_DEFAULT = 64,     &
    LMAX_DEFAULT = 128,      &
    MMAX_DEFAULT = 128

  !---------------------------------------

  ! Compact object parameters.

  ! comb_type
  character(len=S2_STRING_LEN) :: comb_type
  integer, parameter :: COMB_TYPE_OPT_NUM = 7
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TYPE_OPT(COMB_TYPE_OPT_NUM) = (/ &
    'butterfly        ', &
    'gaussian         ', &
    'mexhat           ', &
    'morlet           ', &
    'point            ', &
    'cos_thetaon2     ', &
    'harmonic_gaussian' &
    /)
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TYPE_DEFAULT = COMB_TYPE_OPT(1)

  ! comb_num
  integer :: comb_num
  integer, parameter :: COMB_NUM_OPT(2) = (/ 0, 50 /) 
  integer, parameter :: COMB_NUM_DEFAULT = 10

  ! comb_seed
  integer :: comb_seed

  ! comb_amp_lower
  real(s2_sp) :: comb_amp_lower
  real(s2_sp), parameter :: COMB_AMP_LOWER_OPT(2) = (/ 0.01e0, 5.0e0 /) 
  real(s2_sp), parameter :: COMB_AMP_LOWER_DEFAULT = 0.1e0

  ! comb_amp_upper
  real(s2_sp) :: comb_amp_upper
  real(s2_sp), parameter :: COMB_AMP_UPPER_OPT(2) = (/ 0.01e0, 5.0e0 /) 
  real(s2_sp), parameter :: COMB_AMP_UPPER_DEFAULT = 2.0e0
 
  ! comb_dil_lower
  real(s2_sp) :: comb_dil_lower
  real(s2_sp), parameter :: COMB_DIL_LOWER_OPT(2) = (/ 0.01e0, 1.0e0 /) 
  real(s2_sp), parameter :: COMB_DIL_LOWER_DEFAULT = 0.05e0

  ! comb_dil_upper
  real(s2_sp) :: comb_dil_upper
  real(s2_sp), parameter :: COMB_DIL_UPPER_OPT(2) = (/ 0.01e0, 1.0e0 /) 
  real(s2_sp), parameter :: COMB_DIL_UPPER_DEFAULT = 0.1e0

  ! comb_dil_grid_num
  ! Required so can have exactly same dilations considered in wavelet analysis.
  integer :: comb_dil_grid_num
  integer, parameter :: COMB_DIL_GRID_NUM_DEFAULT = 5

  ! comb_gamma_grid_num
  ! Required so can have exactly same directions considered in wavelet 
  ! analysis.
  integer :: comb_gamma_grid_num
  integer, parameter :: COMB_GAMMA_GRID_NUM_DEFAULT = 5

  ! comb_theta_centre_fov
  real(s2_sp) :: comb_theta_centre_fov
  real(s2_sp), parameter :: COMB_THETA_CENTRE_FOV_OPT(2) = (/ 0.0e0, 180.0e0 /) 
  real(s2_sp), parameter :: COMB_THETA_CENTRE_FOV_DEFAULT = 0e0

  ! comb_tmpl_param_file
  character(len=S2_STRING_LEN) :: comb_tmpl_param_file
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TMPL_FILE_DEFAULT = 'NA'
  real(s2_sp), allocatable :: comb_tmpl_params(:)
  integer :: fileid_tmpl_param = 21
  integer :: n_tmpl_params, iparam
  character(len=1), parameter :: COMMENT_CHAR = '#'
  logical :: params_present = .false.

  !---------------------------------------
  
  ! Primordial CMB parameters.

  ! cmb_status
  logical :: cmb_status
  logical, parameter :: CMB_STATUS_DEFAULT = .true.

  ! cmb_seed
  integer :: cmb_seed

  ! cmb_cl_file
  character(len=S2_STRING_LEN) :: cmb_cl_file
  character(len=S2_STRING_LEN), parameter :: &
    CMB_CL_FILE_DEFAULT = 'data_in/wmap_lcdm_pl_model_yr1_v1.txt'

  !---------------------------------------

  ! Beam parameters

  ! beam_status
  logical :: beam_status
  logical, parameter :: BEAM_STATUS_DEFAULT = .false.

  ! beam_gaussian
  logical :: beam_gaussian
  logical, parameter :: BEAM_GAUSSIAN_DEFAULT = .true.

  ! beam_fwhm
  real(s2_sp) :: beam_fwhm
  real(s2_sp), parameter :: BEAM_FWHM_DEFAULT = 10.0e0

  ! beam_cl_file
  character(len=S2_STRING_LEN) :: beam_cl_file
  character(len=S2_STRING_LEN), parameter :: &
    BEAM_CL_FILE_DEFAULT = 'data_in/map_w1_ampl_bl_yr1_v1.txt'

  !---------------------------------------

  ! Noise parameters

  ! noise_status
  logical :: noise_status
  logical, parameter :: NOISE_STATUS_DEFAULT = .false.

  ! noise_full_sky
  logical :: noise_full_sky
  logical, parameter :: NOISE_FULL_SKY_DEFAULT = .false.

  ! noise_seed
  integer :: noise_seed

  ! noise_std
  real(s2_sp) :: noise_std
  real(s2_sp), parameter :: NOISE_STD_DEFAULT = 0.1e0

  ! noise_sigma0
  real(s2_sp) :: noise_sigma0
  real(s2_sp), parameter :: NOISE_SIGMA0_DEFAULT = 1.0e0

  ! noise_nobs_file
  character(len=S2_STRING_LEN) :: noise_nobs_file
  character(len=S2_STRING_LEN), parameter :: &
    NOISE_NOBS_FILE_DEFAULT = 'data_in/wmap_w1_cleanimap_yr1_v1.fits'

  ! noise_nobs_extension
  integer :: noise_nobs_extension
  integer, parameter :: NOISE_NOBS_EXTENSION_DEFAULT = 2

  !---------------------------------------

  ! Output filename parameters

 ! comb_type
  character(len=S2_STRING_LEN) :: file_type
  integer :: file_type_int
  integer, parameter :: FILE_TYPE_OPT_NUM = 3
  character(len=S2_STRING_LEN), parameter :: &
    FILE_TYPE_OPT(FILE_TYPE_OPT_NUM) = (/ &
    'map', &
    'alm', &
    'sky' &
    /)
  character(len=S2_STRING_LEN), parameter :: &
    FILE_TYPE_DEFAULT = FILE_TYPE_OPT(1)

  ! output_file_csky
  character(len=S2_STRING_LEN) :: output_file_csky
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_CSKY_DEFAULT = 'data_out/csky.fits'

  ! output_file_obj
  character(len=S2_STRING_LEN) :: output_file_obj
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_OBJ_DEFAULT = 'data_out/obj.fits'

  ! output_file_obj_mother
  character(len=S2_STRING_LEN) :: output_file_obj_mother
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_OBJ_MOTHER_DEFAULT = 'data_out/obj_mother.fits'

  ! output_file_cmb
  character(len=S2_STRING_LEN) :: output_file_cmb
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_CMB_DEFAULT = 'data_out/cmb.fits'

  ! output_file_wnoise
  character(len=S2_STRING_LEN) :: output_file_wnoise
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_WNOISE_DEFAULT = 'data_out/wnoise.fits'

  ! output_file_wnoise_nobs
  character(len=S2_STRING_LEN) :: output_file_wnoise_nobs
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_WNOISE_NOBS_DEFAULT = 'data_out/wnoise_nobs.fits'

  ! output_file_wnoise_std
  character(len=S2_STRING_LEN) :: output_file_wnoise_std
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_WNOISE_STD_DEFAULT = 'data_out/wnoise_std.fits'

  ! output_file_param
  character(len=S2_STRING_LEN) :: output_file_param
  character(len=S2_STRING_LEN), parameter :: &
    OUTPUT_FILE_PARAM_DEFAULT = 'data_out/param_out.par'

  !---------------------------------------

  ! Data vairables

  real(s2_sp) :: comb_gamma_lower, comb_gamma_upper
  integer :: iopt, iobj,  fail = 0

  real(s2_sp), allocatable :: amplitude(:)
  real(s2_sp), allocatable :: dilation(:)
  real(s2_sp), allocatable :: alpha(:), beta(:), gamma(:)

  real(s2_sp) :: theta_centre_fov = 0.0
  real(s2_sp) :: xtmp, ytmp, rtmp, phi0, theta0, xmax_sgp

  type(comb_obj) :: obj_mother
  type(comb_csky) :: csky
  type(s2_pl) :: beam
  type(s2_cmb) :: cmb
  type(s2_wnoise) :: wnoise

  !---------------------------------------

  ! Interfaces of external funcitons.

  interface comb_csim_param_check_type

    function comb_csim_param_check_type_char(param, types) result(pass)
      character(len=*), intent(in) :: param
      character(len=*), intent(in) :: types(:)
      logical :: pass
    end function comb_csim_param_check_type_char

    function comb_csim_param_check_type_int(param, range) result(pass)
      integer, intent(in) :: param
      integer, intent(in) :: range(2)
      logical :: pass
    end function comb_csim_param_check_type_int
    
    function comb_csim_param_check_type_real(param, range) result(pass)
      real, intent(in) :: param
      real, intent(in) :: range(2)
      logical :: pass
    end function comb_csim_param_check_type_real
 
  end interface

  interface

     subroutine comb_csim_snap_to_grid(vals, lower, upper, num)
       use s2_types_mod
       real(s2_sp), intent(inout) :: vals(:)
       real(s2_sp), intent(in) :: lower, upper
       integer, intent(in) :: num
     end subroutine comb_csim_snap_to_grid
     
   end interface

  !----------------------------------------------------------------------------

  ! Initialise program.

  call comb_csim_print_header()

  call comb_csim_param_read()

!  call comb_csim_param_write()

  ! Allocate space for parameters.
  allocate(amplitude(comb_num), stat=fail)
  allocate(dilation(comb_num), stat=fail)
  allocate(alpha(comb_num), stat=fail)
  allocate(beta(comb_num), stat=fail)
  allocate(gamma(comb_num), stat=fail)
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_csim')
  end if

  call comb_csim_param_gen()

  ! Centre object if only one.
  if(comb_num == 1) then
     alpha = 0.0e0
     beta = pi/2.0e0
     gamma = 0.0e0
  end if

  ! Read template parameters from file
  if (trim(comb_tmpl_param_file) /= trim(COMB_TMPL_FILE_DEFAULT)) then

     ! Open file
     open(fileid_tmpl_param, file=comb_tmpl_param_file, &
          form='formatted', status='old')
     
     ! Ignore leading comment lines.
     line = COMMENT_CHAR
     do while(line(1:1) == COMMENT_CHAR)
        read(fileid_tmpl_param,'(a)') line
     end do

     ! Read number of source positions from last line read (that is actually
     ! not a comment line).
     read(line, *) line2, n_tmpl_params

     ! Allocate space for parameters.
     allocate(comb_tmpl_params(1:n_tmpl_params), stat=fail)
     if(fail /= 0) then
        call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_csim')
     end if

     ! Read data for sources found.
     do iparam = 1,n_tmpl_params
        read(fileid_tmpl_param,*) comb_tmpl_params(iparam)
     end do

     ! Close file
     close(fileid_tmpl_param)

     ! Set param status.
     params_present = .true.

  end if

  ! ---------------------------------------

  ! Generate simulations.


  !---------------------------------------
  ! Read in beam
  !---------------------------------------

  if(beam_status) then
     if(beam_gaussian) then
        beam = s2_pl_init_guassian(beam_fwhm, lmax)
     else
        beam = s2_pl_init(beam_cl_file, lmin, lmax, &
          NUM_COMMENT_LINES_BEAM_FILE, scale_in=.false.)
     end if
  end if


  !---------------------------------------
  ! Generate cmb
  !---------------------------------------

  if(cmb_status) then
     if(beam_status) then
       cmb = s2_cmb_init(cmb_cl_file, nside, lmin, lmax, &
         NUM_COMMENT_LINES_CLT_FILE, cmb_seed, beam) 
     else
       cmb = s2_cmb_init(cmb_cl_file, nside, lmin, lmax, &
         NUM_COMMENT_LINES_CLT_FILE, cmb_seed) 
     end if
  end if


  !---------------------------------------
  ! Generate white noise
  !---------------------------------------

  if(noise_status) then

    if(noise_full_sky) then
      wnoise = s2_wnoise_init(noise_nobs_file, noise_nobs_extension, &
        noise_sigma0, noise_seed)
      if(s2_wnoise_get_nside(wnoise) /= nside) then
        call s2_wnoise_downsample(wnoise, nside)
      end if
    else
      wnoise = s2_wnoise_init(noise_std, nside, noise_seed)    
    end if

! Noise modelling instrumental effects, thus added after convolution
! with beam, i.e. don't convolve with beam here.
!    if(beam_status) then
!       call s2_wnoise_compute_alm(wnoise, lmax, mmax)
!       call s2_wnoise_conv(wnoise, beam)
!    end if

  end if


  !---------------------------------------
  ! Generate csky
  !---------------------------------------

  ! Generate mother object template.
  select case (comb_type)

     case(COMB_TYPE_OPT(1))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_butterfly, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_butterfly, nside, 1.0e0, &
                name=trim(comb_type))
        end if

     case(COMB_TYPE_OPT(2))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_gaussian, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_gaussian, nside, 1.0e0, &
                name=trim(comb_type))
        end if

     case(COMB_TYPE_OPT(3))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_mexhat, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_mexhat, nside, 1.0e0, &
                name=trim(comb_type))
        end if

     case(COMB_TYPE_OPT(4))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_morlet, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_morlet, nside, 1.0e0, &
                name=trim(comb_type))

        end if

     case(COMB_TYPE_OPT(5))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_point, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_point, nside, 1.0e0, &
                name=trim(comb_type))
        end if

     case(COMB_TYPE_OPT(6))
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplmap_cos_thetaon2, nside, 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplmap_cos_thetaon2, nside, 1.0e0, &
                name=trim(comb_type))
        end if

     case(COMB_TYPE_OPT(7))
!TODO: set parameters
!TODO: to define with different sigmas will need to defined full array rathter than through mother
        if (params_present) then
           obj_mother = comb_obj_init(comb_tmplalm_gaussian, lmax, .false., 1.0e0, &
                name=trim(comb_type), param=comb_tmpl_params)
        else
           obj_mother = comb_obj_init(comb_tmplalm_gaussian, lmax, .false., 1.0e0, &
                name=trim(comb_type))
        end if

  end select

  ! Generate csky with cmb and wnoise if specified.
  if(beam_status) then
     
     if(cmb_status .and. noise_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, cmb, wnoise, beam, lmax, mmax)
        
     else if(cmb_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, cmb, beam=beam, lmax=lmax, mmax=mmax)
        
     else if(noise_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, wnoise=wnoise, beam=beam, lmax=lmax, mmax=mmax)   
        
     else
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, beam=beam, lmax=lmax, mmax=mmax)
        
     end if

  else

     if(cmb_status .and. noise_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, cmb, wnoise)
        
     else if(cmb_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, cmb)
        
     else if(noise_status) then
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma, wnoise=wnoise)   
        
     else
        
        csky = comb_csky_init(obj_mother, comb_num, amplitude, dilation, &
             alpha, beta, gamma)
        
     end if


  end if


  !---------------------------------------
  ! Generate output files
  !---------------------------------------

  select case (file_type)

     case(FILE_TYPE_OPT(1))
        file_type_int = S2_SKY_FILE_TYPE_MAP

     case(FILE_TYPE_OPT(2))
        file_type_int = S2_SKY_FILE_TYPE_ALM

     case(FILE_TYPE_OPT(3))
        file_type_int = S2_SKY_FILE_TYPE_SKY

  end select

  call comb_csky_write_sky_full(csky, output_file_csky, file_type_in=file_type_int)
  call comb_csky_write_sky_obj(csky, output_file_obj, file_type_in=file_type_int)
  call comb_obj_write_sky(obj_mother, output_file_obj_mother, file_type_in=file_type_int)

  if(cmb_status) then
     call s2_cmb_write_sky(cmb, trim(output_file_cmb))
  end if

  if(noise_status) then
     call s2_wnoise_write_sky_file(wnoise, output_file_wnoise)
     if(noise_full_sky) then
      call s2_wnoise_write_nobs_file(wnoise, output_file_wnoise_nobs)
      call s2_wnoise_write_std_file(wnoise, output_file_wnoise_std)
     end if
  end if

  call comb_csky_write_param(csky, output_file_param)


  !---------------------------------------
  ! Free memory
  !---------------------------------------

  if(comb_csky_get_init(csky)) call comb_csky_free(csky)
  if(comb_obj_get_init(obj_mother)) call comb_obj_free(obj_mother)
  if(s2_pl_get_init(beam)) call s2_pl_free(beam)
  if(s2_cmb_get_init(cmb)) call s2_cmb_free(cmb)
  if(s2_wnoise_get_init(wnoise)) call s2_wnoise_free(wnoise)
  deallocate(amplitude)
  deallocate(dilation)
  deallocate(alpha, beta, gamma)


  !----------------------------------------------------------------------------

  contains 


    !--------------------------------------------------------------------------
    ! comb_csim_print_header
    ! 
    !! Print comb_csim header to standard out.
    !!
    !! Variables:
    !!   - Has local scope so accesses variables defined in comb_csim program.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_csim_print_header

      write(*,'(a,a)') '***************************************************', &
        '********************'
      write(*,'(a,a,a)') '* ', CODE_NAME, &
        '                                                           *'
      write(*,'(a,a,a)') '* Version ', CODE_VERSION, &
        '                                                         *'
      write(*,'(a,a,a)') '* ', CODE_DESCRIPTION, &
        '                              *'
      write(*,'(a,a,a)') '* Written by ', CODE_AUTHOR, &
        '                                             *'
      write(*,'(a,a)') '***************************************************', &
           '********************'

    end subroutine comb_csim_print_header


    !--------------------------------------------------------------------------
    ! comb_csim_param_read
    !    
    !! Parse comb_csim input parameters.  Either read from input parameter
    !! file or read interactively from standard input.
    !!
    !! Variables:
    !!   - Has local scope so accesses variables defined in comb_csim program.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_csim_param_read

      ! Initialise file parser.
      if(nArguments() == 0) then
         filename_param = ''
      else
         if(nArguments() /= 1) then
            call comb_error(COMB_ERROR_CSIM_NARG, 'comb_csim', &
            comment_add='Usage: comb_csim [input parameter filename]')
         end if
         call getArgument(1, filename_param)
      end if
      handle = parse_init(trim(filename_param))
  
      !---------------------------------------

      ! Get Healpix parameters.

      ! Get nside.
      description = concatnl("", &
        "Enter the resolution parameter (nside) for the simulated sky: ", &
        "(npix = 12*nside**2, where nside must be a power of 2)")
3     continue
      nside = parse_int(handle, 'nside', &
        default=NSIDE_DEFAULT, descr=description)
      if(nside2npix(nside) < 0) then
         if(handle%interactive) goto 3
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_csim', &
              comment_add='nside invalid')
      endif

      ! Get lmax.  
      write(line,"(a,i5,a)") "Require: (0 < l_max <= ",3*nside-1,")"
      description = concatnl(&
        & "", &
        & "Enter the harmonic resolution lmax: ", &
        & line )
4     continue
      lmax = parse_int(handle, 'lmax', &
        default=LMAX_DEFAULT, descr=description)
      if(lmax <= 0 .or. lmax >= 3*nside) then
         if(handle%interactive) goto 4
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_csim', &
           comment_add='lmax invalid')
      end if

      ! Get mmax.  
      write(line,"(a,i5,a)") "Require: (0 < mmax <= lmax = ", lmax, ")"
      description = concatnl(&
        & "", &
        & "Enter the harmonic resolution mmax: ", &
        & line)
5     continue
      mmax = parse_int(handle, 'mmax', &
        default=MMAX_DEFAULT, descr=description)
      if(mmax <= 0 .or. mmax >= 3*nside) then
         if(handle%interactive) goto 5
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_csim', &
           comment_add='lmax invalid')
      end if

      !---------------------------------------

      ! Get compact object (cpac) parameters.

      ! Get comb_type.
      line = '(Options:'
      iopt = 1
      line = trim(line) // ' ' // trim(COMB_TYPE_OPT(iopt)) // ';'
      do iopt = 2,COMB_TYPE_OPT_NUM
         line = trim(line) // ' ' // trim(COMB_TYPE_OPT(iopt))
         if(iopt /= COMB_TYPE_OPT_NUM) line = trim(line) // ';'
      end do
      line = trim(line) // ')'
      description = concatnl('', &
        'Enter comb_type: ', &
        line)
6     continue
      comb_type = parse_string(handle, 'comb_type', &
        default=trim(COMB_TYPE_DEFAULT), descr=description)
      if(.not. comb_csim_param_check_type(comb_type, COMB_TYPE_OPT)) then
         if(handle%interactive) goto 6
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_type invalid')
      end if
      
      ! Get comb_num.
      write(line,'(a,i5,a,i5,a)') '(Between the range ', &
        COMB_NUM_OPT(1), ' and ', COMB_NUM_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_num: ', &
        line)
7     continue
      comb_num = parse_int(handle, 'comb_num', &
        default=COMB_NUM_DEFAULT, descr=description)
      if(.not. comb_csim_param_check_type(comb_num, COMB_NUM_OPT)) then
         if(handle%interactive) goto 7
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_num invalid')
      end if

      ! Get comb_seed.
      description = concatnl('', &
        'Enter comb_seed: ')
      call system_clock(comb_seed)
      comb_seed = parse_int(handle, 'comb_seed', &
        default=comb_seed, descr=description)

      ! Get comb_amp_lower.
      write(line,'(a,f4.1,a,f4.1,a)') '(Between the range ', &
        COMB_AMP_LOWER_OPT(1), ' and ', COMB_AMP_LOWER_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_amp_lower: ', &
        line)
8     continue
      comb_amp_lower = parse_real(handle, 'comb_amp_lower', &
        default=COMB_AMP_LOWER_DEFAULT, descr=description)
      if(.not. &
        comb_csim_param_check_type(comb_amp_lower, COMB_AMP_LOWER_OPT)) then
         if(handle%interactive) goto 8
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_amp_lower invalid')
      end if

      ! Get comb_amp_upper.
      write(line,'(a,f4.1,a,f4.1,a)') '(Between the range ', &
        COMB_AMP_UPPER_OPT(1), ' and ', COMB_AMP_UPPER_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_amp_upper: ', &
        line)
9     continue
      comb_amp_upper = parse_real(handle, 'comb_amp_upper', &
        default=COMB_AMP_UPPER_DEFAULT, descr=description)
      if(.not. &
        comb_csim_param_check_type(comb_amp_upper, COMB_AMP_UPPER_OPT) &
        .or. comb_amp_upper < comb_amp_lower) then
         if(handle%interactive) goto 9
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_amp_upper invalid')
      end if

      ! Get comb_dil_lower.
      write(line,'(a,f4.1,a,f4.1,a)') '(Between the range ', &
        COMB_DIL_LOWER_OPT(1), ' and ', COMB_DIL_LOWER_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_dil_lower: ', &
        line)
10    continue
      comb_dil_lower = parse_real(handle, 'comb_dil_lower', &
        default=COMB_DIL_LOWER_DEFAULT, descr=description)
      if(.not. &
        comb_csim_param_check_type(comb_dil_lower, COMB_DIL_LOWER_OPT)) then
         if(handle%interactive) goto 10
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_dil_lower invalid')
      end if

      ! Get comb_dil_upper.
      write(line,'(a,f4.1,a,f4.1,a)') '(Between the range ', &
        COMB_DIL_UPPER_OPT(1), ' and ', COMB_DIL_UPPER_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_dil_upper: ', &
        line)
11    continue
      comb_dil_upper = parse_real(handle, 'comb_dil_upper', &
        default=COMB_DIL_UPPER_DEFAULT, descr=description)
      if(.not. &
        comb_csim_param_check_type(comb_dil_upper, COMB_DIL_UPPER_OPT) &
        .or. comb_dil_upper < comb_dil_lower) then
         if(handle%interactive) goto 11
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_dil_upper invalid')
      end if

      ! Get comb_dil_grid_num.
      write(line,'(a,a,a)') &
        '(In order to construct the allowable equispaced dilation grid. ', &
        'Sample dilations are randomly selected from the allowable grid. ', &
        'Set to -1 for continuous         dilation range.)'
      description = concatnl('', &
        'Enter comb_dil_grid_num: ', &
        line)
12    continue
      comb_dil_grid_num = parse_int(handle, 'comb_dil_grid_num', &
        default=COMB_DIL_GRID_NUM_DEFAULT, descr=description)
      if( comb_dil_grid_num < -1 .or. comb_dil_grid_num == 0 &
           .or. (comb_dil_grid_num == 1 &
                 .and. abs(comb_dil_lower - comb_dil_upper) > TOL) ) then
         if(handle%interactive) goto 12
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_dil_grid_num invalid')
      end if

      ! Get comb_gamma_grid_num.
      write(line,'(a,a,a)') &
        '(In order to construct the allowable equispaced gamma grid. ', &
        'Sample angles are randomly selected from the allowable grid. ', &
        'Set to -1 for continuous         gamma range.)'
      description = concatnl('', &
        'Enter comb_gamma_grid_num: ', &
        line)
13    continue
      comb_gamma_grid_num = parse_int(handle, 'comb_gamma_grid_num', &
        default=COMB_GAMMA_GRID_NUM_DEFAULT, descr=description)
      if( comb_gamma_grid_num < -1 .or. comb_gamma_grid_num == 0 ) then !&
!           .or.comb_gamma_grid_num == 1 ) then
         if(handle%interactive) goto 13
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_gamma_grid_num invalid')
      end if

      ! Get comb_theta_centre_fov.
      write(line,'(a,f4.1,a,f4.1,a)') '(Between the range ', &
        COMB_THETA_CENTRE_FOV_OPT(1), ' and ', COMB_THETA_CENTRE_FOV_OPT(2), ')'
      description = concatnl('', &
        'Enter comb_theta_centre_fov (in degrees): ', &
        line)
131   continue
      comb_theta_centre_fov = parse_real(handle, 'comb_theta_centre_fov', &
        default=COMB_THETA_CENTRE_FOV_DEFAULT, descr=description)
      if(.not. &
        comb_csim_param_check_type(comb_theta_centre_fov, COMB_THETA_CENTRE_FOV_OPT)) then
         if(handle%interactive) goto 131
         call comb_error(COMB_ERROR_CSIM_PARAM_INVALID, 'comb_sim', &
              comment_add='comb_theta_centre_fov invalid')
      end if
      ! Convert to radians.
      comb_theta_centre_fov = comb_theta_centre_fov / 180.0 * pi

      ! Get comb_tmpl_param_file.
      description = concatnl('', &
        'Enter comb_tmpl_param_file (set to "NA" if template parameter overrides not required): ')
      comb_tmpl_param_file = parse_string(handle, 'comb_tmpl_param_file', &
        default=trim(COMB_TMPL_FILE_DEFAULT), descr=description)

      !---------------------------------------

      ! Get cmb parameters.

      ! Get cmb_status
      description = concatnl(&
        & "", &
        & "Enter cmb_status: ")
      cmb_status = parse_lgt(handle, 'cmb_status', &
        default=CMB_STATUS_DEFAULT, descr=description)

      ! Get comb_seed.
      description = concatnl('', &
        'Enter cmb_seed: ')
      call system_clock(cmb_seed)
      cmb_seed = parse_int(handle, 'cmb_seed', &
        default=cmb_seed, descr=description)

      ! Get cmb_cl_file.
      description = concatnl('', &
        'Enter cmb_cl_file: ')
      cmb_cl_file = parse_string(handle, 'cmb_cl_file', &
        default=trim(CMB_CL_FILE_DEFAULT), descr=description)

      !---------------------------------------

      ! Get beam parameters.

      ! Get beam_status.
      description = concatnl(&
        & "", &
        & "Enter beam_status: ")
      beam_status = parse_lgt(handle, 'beam_status', &
        default=BEAM_STATUS_DEFAULT, descr=description)

      ! Get beam_gaussian.
      description = concatnl(&
        & "", &
        & "Enter beam_gaussian: ")
      beam_gaussian = parse_lgt(handle, 'beam_gaussian', &
        default=BEAM_GAUSSIAN_DEFAULT, descr=description)

      if(beam_gaussian) then
         ! Get beam fwhm.
         description = concatnl('', &
           'Enter beam_fwhm (arcmin): ')
         beam_fwhm = parse_real(handle, 'beam_fwhm', &
           default=BEAM_FWHM_DEFAULT, descr=description)
         beam_fwhm = s2_vect_arcmin_to_rad(beam_fwhm)
      else
         ! Get beam_cl_file.
         description = concatnl('', &
           'Enter beam_cl_file: ')
         beam_cl_file = parse_string(handle, 'beam_cl_file', &
           default=trim(BEAM_CL_FILE_DEFAULT), descr=description)
      end if

      !---------------------------------------

      ! Get noise parameters.

      ! Get noise_status.
      description = concatnl(&
        & "", &
        & "Enter noise_status: ")
      noise_status = parse_lgt(handle, 'noise_status', &
        default=NOISE_STATUS_DEFAULT, descr=description)

      ! Get noise_full_sky.
      description = concatnl(&
        & "", &
        & "Enter noise_full_sky: ")
      noise_full_sky = parse_lgt(handle, 'noise_full_sky', &
        default=NOISE_FULL_SKY_DEFAULT, descr=description)

      ! Get noise_seed.
      description = concatnl('', &
        'Enter noise_seed: ')
      call system_clock(noise_seed)
      noise_seed = parse_int(handle, 'noise_seed', &
        default=noise_seed, descr=description)

      ! Get noise_std.
      description = concatnl('', &
        'Enter noise_std: ')
      noise_std = parse_real(handle, 'noise_std', &
        default=NOISE_STD_DEFAULT, descr=description)

      ! Get noise_sigma0.
      description = concatnl('', &
        'Enter noise_sigma0: ')
      noise_sigma0 = parse_real(handle, 'noise_sigma0', &
        default=NOISE_SIGMA0_DEFAULT, descr=description)

      ! Get noise_nobs_file.
      description = concatnl('', &
        'Enter noise_nobs_file: ')
      noise_nobs_file = parse_string(handle, 'noise_nobs_file', &
        default=trim(NOISE_NOBS_FILE_DEFAULT), descr=description)

      ! Get noise_nobs_extension.
      description = concatnl('', &
        'Enter noise_nobs_extension: ')
      noise_nobs_extension = parse_int(handle, 'noise_nobs_extension', &
        default=NOISE_NOBS_EXTENSION_DEFAULT, descr=description)

      !---------------------------------------

      ! Get output filename parameters.

      ! Get file_type.
      line = '(Options:'
      iopt = 1
      line = trim(line) // ' ' // trim(FILE_TYPE_OPT(iopt)) // ';'
      do iopt = 2,FILE_TYPE_OPT_NUM
         line = trim(line) // ' ' // trim(FILE_TYPE_OPT(iopt))
         if(iopt /= FILE_TYPE_OPT_NUM) line = trim(line) // ';'
      end do
      line = trim(line) // ')'
      description = concatnl('', &
        'Enter file_type: ', &
        line)
14    continue
      file_type = parse_string(handle, 'file_type', &
        default=trim(FILE_TYPE_DEFAULT), descr=description)
      if(.not. comb_csim_param_check_type(file_type, FILE_TYPE_OPT)) then
         if(handle%interactive) goto 14
         call s2_error(S2_ERROR_SKY_FTYPE_INVALID, 'comb_sim', &
              comment_add='file_type invalid')
      end if

      ! Get ouput_file_csky.
      description = concatnl('', &
        'Enter output_file_csky: ')
      output_file_csky = parse_string(handle, 'output_file_csky', &
        default=trim(OUTPUT_FILE_CSKY_DEFAULT), descr=description)

      ! Get ouput_file_obj.
      description = concatnl('', &
        'Enter output_file_obj: ')
      output_file_obj = parse_string(handle, 'output_file_obj', &
        default=trim(OUTPUT_FILE_OBJ_DEFAULT), descr=description)

      ! Get ouput_file_obj_mother.
      description = concatnl('', &
        'Enter output_file_obj_mother: ')
      output_file_obj_mother = parse_string(handle, 'output_file_obj_mother', &
        default=trim(OUTPUT_FILE_OBJ_MOTHER_DEFAULT), descr=description)

      ! Get ouput_file_cmb.
      description = concatnl('', &
        'Enter output_file_cmb: ')
      output_file_cmb = parse_string(handle, 'output_file_cmb', &
        default=trim(OUTPUT_FILE_CMB_DEFAULT), descr=description)

      ! Get ouput_file_wnoise.
      description = concatnl('', &
        'Enter output_file_wnoise: ')
      output_file_wnoise = parse_string(handle, 'output_file_wnoise', &
        default=trim(OUTPUT_FILE_WNOISE_DEFAULT), descr=description)

      ! Get ouput_file_wnoise_nobs.
      description = concatnl('', &
        'Enter output_file_wnoise_nobs: ')
      output_file_wnoise_nobs = parse_string(handle, &
        'output_file_wnoise_nobs', &
        default=trim(OUTPUT_FILE_WNOISE_NOBS_DEFAULT), descr=description)

      ! Get ouput_file_wnoise_std.
      description = concatnl('', &
        'Enter output_file_wnoise_std: ')
      output_file_wnoise_std = parse_string(handle, 'output_file_wnoise_std', &
        default=trim(OUTPUT_FILE_WNOISE_STD_DEFAULT), descr=description)
     
      ! Get ouput_file_param.
      description = concatnl('', &
        'Enter output_file_param: ')
      output_file_param = parse_string(handle, 'output_file_param', &
        default=trim(OUTPUT_FILE_PARAM_DEFAULT), descr=description)

      write(*,*)

    end subroutine comb_csim_param_read


    !--------------------------------------------------------------------------
    ! comb_csim_param_write
    !
    !! Write parameter values to screen.
    !!
    !! Variables:
    !!   - Has local scope so accesses variables defined in comb_csim program.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_csim_param_write

      write(*,'(a)') 'Writing parameters:'

      write(*,'(a,i10)') ' nside          = ', nside
      write(*,'(a,i10)') ' lmax           = ', lmax
      write(*,'(a,i10)') ' mmax           = ', mmax

      write(*,'(a,i10)') ' comb_seed      = ', comb_seed
      write(*,'(a,a10)') ' comb_type      = ', trim(comb_type)
      write(*,'(a,i10)') ' comb_num       = ', comb_num
      write(*,'(a,f10.1)') ' comb_amp_lower = ', comb_amp_lower
      write(*,'(a,f10.1)') ' comb_amp_upper = ', comb_amp_upper
      write(*,'(a,f10.1)') ' comb_dil_lower = ', comb_dil_lower
      write(*,'(a,f10.1)') ' comb_dil_upper = ', comb_dil_upper
      write(*,'(a,i5)') ' comb_dil_grid_num   = ', comb_dil_grid_num
      write(*,'(a,i5)') ' comb_gamma_grid_num = ', comb_gamma_grid_num

      write(*,'(a,l10)') ' cmb_status     = ', cmb_status
      write(*,'(a,i10)') ' cmb_seed       = ', cmb_seed
      write(*,'(a,a)') ' cmb_cl_file    = ', trim(cmb_cl_file)

      write(*,'(a,l10)') ' beam_status    = ', beam_status
      write(*,'(a,a)') ' beam_cl_file   = ', trim(beam_cl_file)

      write(*,'(a,l10)') ' noise_status   = ', noise_status
      write(*,'(a,l10)') ' noise_full_sky = ', noise_full_sky
      write(*,'(a,i10)') ' noise_seed     = ', noise_seed
      write(*,'(a,f10.1)') ' noise_std      = ', noise_std
      write(*,'(a,f10.1)') ' noise_sigma0   = ', noise_sigma0
      write(*,'(a,a)') ' noise_nobs_file = ', trim(noise_nobs_file)
      write(*,'(a,i4)') ' noise_nobs_extension = ', noise_nobs_extension

      write(*,'(a,a)') ' output_file_csky = ', trim(output_file_csky)
      write(*,'(a,a)') ' output_file_obj = ', trim(output_file_obj)
      write(*,'(a,a)') ' output_file_obj_mother= ',trim(output_file_obj_mother)
      write(*,'(a,a)') ' output_file_cmb = ', trim(output_file_cmb)
      write(*,'(a,a)') ' output_file_wnoise = ', trim(output_file_wnoise)
      write(*,'(a,a)') ' output_file_wnoise_nobs = ', trim(output_file_wnoise_nobs)
      write(*,'(a,a)') ' output_file_wnoise_std = ', trim(output_file_wnoise_std)
      write(*,'(a,a)') ' output_file_param = ', trim(output_file_param)

    end subroutine comb_csim_param_write


    !--------------------------------------------------------------------------
    ! comb_csim_param_gen
    ! 
    !! Generate samples compact object parameters in specified range.
    !!
    !! Notes:
    !!   - For now gamma angle and dilation are snapped to an equispaced grid
    !!     so that corresponds to angles considered in a usual spherical 
    !!     wavelet analysis.
    !!
    !! Variables:
    !!   - Has local scope so accesses variables defined in comb_csim program.
    !
    !! @author J. D. McEwen
    !! @version 0.1 September 2004
    !
    ! Revisions:
    !   September 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine comb_csim_param_gen()

      ! Define gamma range.
      comb_gamma_lower = 0.0e0
      comb_gamma_upper = 2.0e0*pi - 2.0e0*pi/real(comb_gamma_grid_num, s2_sp)

      ! Generate random parameters satisfying uniform distribution in
      ! specified range.
      do iobj = 1,comb_num
     
         amplitude(iobj) = s2_distn_sample_uniform(comb_seed, &
              comb_amp_lower, comb_amp_upper)
         
         dilation(iobj) = s2_distn_sample_uniform(comb_seed, comb_dil_lower, &
              comb_dil_upper)

         if(abs(comb_theta_centre_fov) < TOL) then
            alpha(iobj) = s2_distn_sample_uniform(comb_seed, 0.0e0, 2*pi)
            beta(iobj) = s2_distn_sample_uniform(comb_seed, 0.0e0, pi)
         else
            xmax_sgp = sqrt(2.0) * tan(comb_theta_centre_fov/4.0)
            xtmp = s2_distn_sample_uniform(comb_seed, -xmax_sgp, xmax_sgp)
            ytmp = s2_distn_sample_uniform(comb_seed, -xmax_sgp, xmax_sgp)
            rtmp = sqrt(xtmp**2 + ytmp**2)
            phi0 = atan2(ytmp, xtmp)
            theta0 = 2 * atan(rtmp/2.0)
            alpha(iobj) = mod(phi0, 2*pi)
            beta(iobj) = mod(theta0, pi)
         end if

         gamma(iobj) = s2_distn_sample_uniform(comb_seed, comb_gamma_lower, &
              comb_gamma_upper)
         
      end do
  
      ! Snap dilations and gamma angles to quantised grids that correspond 
      ! to discretised wavelet transform domains (if required).
      if(comb_dil_grid_num /= -1) then
         call comb_csim_snap_to_grid(dilation, comb_dil_lower, &
              comb_dil_upper, comb_dil_grid_num)
      end if
      if(comb_gamma_grid_num /= -1) then
         call comb_csim_snap_to_grid(gamma, comb_gamma_lower, &
              comb_gamma_upper, comb_gamma_grid_num)
      end if
      
    end subroutine comb_csim_param_gen


end program comb_csim


!------------------------------------------------------------------------------
! External routines
!------------------------------------------------------------------------------

!--------------------------------------------------------------------------
! comb_csim_param_check_type_char
!
!! Check string parameter matches one of the types specified in the types
!! string array.
!!
!! Variables:
!!   - param: String parameter.
!!   - types: Array of strings containing allowable types.
!!   - pass: Logical indicating whether the string parameter param matches
!!     one of the strings in the types array.
!
!! @author J. D. McEwen
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!--------------------------------------------------------------------------

function comb_csim_param_check_type_char(param, types) result(pass)
  
  implicit none
  
  character(len=*), intent(in) :: param
  character(len=*), intent(in) :: types(:)
  logical :: pass
  
  integer :: ntypes, i
  
  ntypes = size(types)
  pass = .false.
  
  do i = 1,ntypes
     if(param == types(i)) pass = .true.
  end do
  
end function comb_csim_param_check_type_char


!--------------------------------------------------------------------------
! comb_csim_param_check_type_int
!
!! Check integer parameter lies within the range specified by range.
!!
!! Variables:
!!   - param: Integer parameter.
!!   - range: Integer range of length 2 containing the lower and upper
!!     allowable bounds for the integer parameter.
!!   - pass: Logical indicating whether tht integer parameter lies between
!!     the two specified bounds.
!
!! @author J. D. McEwen
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!--------------------------------------------------------------------------

function comb_csim_param_check_type_int(param, range) result(pass)
  
  implicit none
  
  integer, intent(in) :: param
  integer, intent(in) :: range(2)
  logical :: pass
  
  if(param < range(1) .or. param > range(2)) then
     pass = .false.
  else
     pass = .true.
  end if
  
end function comb_csim_param_check_type_int


!--------------------------------------------------------------------------
! comb_csim_param_check_type_real
! 
!! Check real parameter lies within the range specified by range.
!!
!! Variables:
!!   - param: Real parameter.
!!   - range: Real range of length 2 containing the lower and upper
!!     allowable bounds for the real parameter.
!!   - pass: Logical indicating whether tht real parameter lies between the
!!     two specified bounds.
!
!! @author J. D. McEwen
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!--------------------------------------------------------------------------

function comb_csim_param_check_type_real(param, range) result(pass)
  
  implicit none
  
  real, intent(in) :: param
  real, intent(in) :: range(2)
  logical :: pass
  
  if(param < range(1) .or. param > range(2)) then
     pass = .false.
  else
     pass = .true.
  end if
  
end function comb_csim_param_check_type_real


!--------------------------------------------------------------------------
! comb_csim_snap_to_grid
!
!! Snap value array (vals) to posisions that lie on an equisampled grid. 
!! Required so that have parameter values that correspond to positions
!! directly probed with usual wavelet analysis.
!!
!! Notes:
!!   - All grid positions should be equally likely.  This inc2 specifies
!!     the space between sample positions and inc1 specifies spacing for
!!     snapping values to the grid.  See diagram (for num_grid=5):
!!
!!     l.................. u 
!!
!!     [...[...[...[...[......(spacings to snap positions to grid points)
!!
!!     |----|----|----|----|..(sample positions)
!!
!!     <--->...(inc1)
!!
!!     <---->..(inc2)
!!
!! Variables:
!!   - vals:
!!   - lower: Lower bound of values.
!!   - upper: upper bound of values.
!!   - num_grid: Number of grid samples between lower and upper bounds.
!
!! @author J. D. McEwen
!! @version 0.1 September 2004
!
! Revisions:
!   September 2004 - Written by Jason McEwen
!--------------------------------------------------------------------------

subroutine comb_csim_snap_to_grid(vals, lower, upper, num_grid)
  
  use s2_types_mod
  
  implicit none
  
  real(s2_sp), intent(inout) :: vals(:)
  real(s2_sp), intent(in) :: lower, upper
  integer :: num_grid
  
  integer :: i, n
  real(s2_sp) :: inc1, inc2
  real(s2_sp) :: TOL = 1e-5


  ! Do nothing if lower equals upper.
  if(abs(lower-upper) < TOL) return

  n = size(vals)
  inc1 = (upper - lower) / real(num_grid-1, s2_sp)
  inc2 = (upper - lower) / real(num_grid, s2_sp)

  do i = 1,n
     vals(i) = floor((vals(i) - lower) / inc2) * inc1 + lower
  end do
  
end subroutine comb_csim_snap_to_grid

