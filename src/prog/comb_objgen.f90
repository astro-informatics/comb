!------------------------------------------------------------------------------
! comb_objgen -- COMB program
!
!! Product a sky map containing a number of template objects with various 
!! parameters.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 May 2006
!
! Revisions:
!   May 2006 - Written by Jason McEwen
!------------------------------------------------------------------------------

program comb_objgen

  use s2_types_mod, only: S2_STRING_LEN, s2_sp, pi
  use s2_error_mod
  use s2_pl_mod
  use comb_error_mod
  use comb_tmplmap_mod
  use comb_obj_mod
  use comb_csky_mod
  
  use extension, only: getArgument, nArguments

  implicit none

  !----------------------------------------------------------------------------

  character(len=1), parameter :: COMMENT_CHAR = '#'
  character(len=S2_STRING_LEN) :: filename_inp
  character(len=S2_STRING_LEN) :: filename_out
  character(len=S2_STRING_LEN) :: line, line2

  integer :: nside = 64, lmax = 128, mmax = 128, fail, fileid = 32

  character(len=S2_STRING_LEN) :: comb_type
  integer, parameter :: COMB_TYPE_OPT_NUM = 6
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TYPE_OPT(COMB_TYPE_OPT_NUM) = (/ &
    'butterfly   ', &
    'gaussian    ', &
    'mexhat      ', &
    'morlet      ', &
    'point       ', &
    'cos_thetaon2'  &
    /)
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TYPE_DEFAULT = COMB_TYPE_OPT(1)

  logical, parameter :: BEAM_STATUS_DEFAULT = .false.
  logical :: beam_status = BEAM_STATUS_DEFAULT

  real(s2_sp), parameter :: BEAM_FWHM_DEFAULT = 13.0e0
  real(s2_sp) :: beam_fwhm = BEAM_FWHM_DEFAULT

  integer :: n_source, i_source
  real(s2_sp) :: dilation = 0.1
  real(s2_sp), allocatable :: dilation_array(:)
  real(s2_sp), allocatable :: amplitude(:)
  real(s2_sp), allocatable :: alpha(:), beta(:), gamma(:)

  type(comb_csky) :: csky
  type(comb_obj) :: obj_mother
  type(s2_pl) :: beam


  !---------------------------------------
  ! Parse input parameters
  !---------------------------------------

  call parse_options()


  !---------------------------------------
  ! Read localised regions from file
  !---------------------------------------

  open(fileid, file=filename_inp, form='formatted', status='old')

  ! Ignore leading comment lines.
  line = COMMENT_CHAR
  do while(line(1:1) == COMMENT_CHAR)
     read(fileid,'(a)') line
  end do

  ! Read number of source positions from last line read (that is actually
  ! not a comment line).
  read(line, *) line2, n_source

  ! Allocate space for parameters.
  allocate(amplitude(1:n_source), stat=fail)
  allocate(alpha(1:n_source), stat=fail)
  allocate(beta(1:n_source), stat=fail)
  allocate(gamma(1:n_source), stat=fail)
  allocate(dilation_array(1:n_source), stat=fail)
  dilation_array(1:n_source) = dilation
  if(fail /= 0) then
     call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_objgen')
  end if

  ! Read data for sources found.
  do i_source = 1,n_source
     read(fileid,'(a)') line   ! Ignore comment line.
     read(fileid,*) line, amplitude(i_source)
     read(fileid,*) line, alpha(i_source)
     read(fileid,*) line, beta(i_source)
     read(fileid,*) line, gamma(i_source)
  end do


  !---------------------------------------
  ! Read in beam
  !---------------------------------------

  if(beam_status) then
     beam = s2_pl_init_guassian(beam_fwhm, lmax)
  end if


  !---------------------------------------
  ! Generate csky
  !---------------------------------------

  ! Initialise mother object.
  select case (comb_type)

     case(COMB_TYPE_OPT(1))
        obj_mother = comb_obj_init(comb_tmplmap_butterfly, nside, 1.0e0, &
             name=trim(comb_type))

     case(COMB_TYPE_OPT(2))
        obj_mother = comb_obj_init(comb_tmplmap_gaussian, nside, 1.0e0, &
             name=trim(comb_type))

     case(COMB_TYPE_OPT(3))
        obj_mother = comb_obj_init(comb_tmplmap_mexhat, nside, 1.0e0, &
             name=trim(comb_type))

     case(COMB_TYPE_OPT(4))
        obj_mother = comb_obj_init(comb_tmplmap_morlet, nside, 1.0e0, &
             name=trim(comb_type))

     case(COMB_TYPE_OPT(5))
        obj_mother = comb_obj_init(comb_tmplmap_point, nside, 1.0e0, &
             name=trim(comb_type))

     case(COMB_TYPE_OPT(6))
        obj_mother = comb_obj_init(comb_tmplmap_cos_thetaon2, nside, 1.0e0, &
             name=trim(comb_type))

  end select

  ! Generate csky.
  if(beam_status) then
     csky = comb_csky_init(obj_mother, n_source, amplitude, &
          dilation_array, alpha, beta, gamma, beam=beam)
  else
     csky = comb_csky_init(obj_mother, n_source, amplitude, &
          dilation_array, alpha, beta, gamma)
  end if


  !---------------------------------------
  ! Generate output files
  !---------------------------------------

  call comb_csky_write_sky_obj(csky, filename_out)


  !---------------------------------------
  ! Free memory
  !---------------------------------------

  call comb_csky_free(csky)
  call comb_obj_free(obj_mother)
  if(beam_status) call s2_pl_free(beam)
  deallocate(amplitude)
  deallocate(dilation_array)
  deallocate(alpha, beta, gamma)


  !----------------------------------------------------------------------------

  contains 


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: n, i
      character(len=S2_STRING_LEN) :: opt
      character(len=S2_STRING_LEN) :: arg
      
      n = nArguments()
     
      do i = 1,n,2
        
        call getArgument(i,opt)
     
        if (i == n .and. trim(opt) /= '-help') then
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: comb_objgen [-inp filename_inp]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-tmpl comb_type]'
            write(*,'(a)') '                   [-dil dilation]'
            write(*,'(a)') '                   [-nside nside]'
            write(*,'(a)') '                   [-lmax lmax (only if beam present)]'
            write(*,'(a)') '                   [-beam beam_fwhm (optional)]'
            stop
          
          case ('-inp')
            filename_inp = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-tmpl')
            read(arg,*) comb_type

          case ('-dil')
            read(arg,*) dilation

          case ('-nside')
            read(arg,*) nside

          case ('-lmax')
            read(arg,*) lmax
            mmax = lmax

          case ('-beam')
            read(arg,*) beam_fwhm
            beam_status = .true.

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program comb_objgen

