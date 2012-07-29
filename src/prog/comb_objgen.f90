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
  integer, parameter :: COMB_TYPE_OPT_NUM = 8
  character(len=S2_STRING_LEN), parameter :: &
    COMB_TYPE_OPT(COMB_TYPE_OPT_NUM) = (/ &
    'butterfly   ', &
    'gaussian    ', &
    'mexhat      ', &
    'morlet      ', &
    'point       ', &
    'cos_thetaon2',  &
    'bubble      ', &
    'texture     ' &
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
  real(s2_sp), allocatable :: source_size(:)

  type(comb_csky) :: csky
  type(comb_obj) :: obj_mother
  type(comb_obj), allocatable :: obj(:)
  type(s2_pl) :: beam

  logical :: include_size = .false.
  real(s2_sp) :: bubble_params(1:4)
  real(s2_sp) :: texture_params(1:1)

  logical :: plot_all = .false.


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
  allocate(source_size(1:n_source), stat=fail)
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
     if(include_size) then
        read(fileid,*) line, source_size(i_source)
        read(fileid,*) line ! Skip significance level.
        read(fileid,*) line ! Skip significance radius.
        read(fileid,*) line ! Skip mask overlap.
     end if
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

     case(COMB_TYPE_OPT(7))
        allocate(obj(1:n_source), stat=fail)
        if(fail /= 0) then
           call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_objgen')
        end if
        do i_source = 1,n_source
           bubble_params(1) = 1.0
           bubble_params(2) = 0.2
           bubble_params(3) = source_size(i_source) / PI * 180
           bubble_params(4) = bubble_params(3) * 1.1
           obj(i_source) = comb_obj_init(comb_tmplmap_bubble, nside, &
                amplitude(i_source), dilation=1.0, &
                alpha=alpha(i_source), beta=beta(i_source), gamma=gamma(i_source), &
                name=trim(comb_type), param=bubble_params)
        end do

     case(COMB_TYPE_OPT(8))
        allocate(obj(1:n_source), stat=fail)
        if(fail /= 0) then
           call s2_error(S2_ERROR_MEM_ALLOC_FAIL, 'comb_objgen')
        end if
        do i_source = 1,n_source
           texture_params(1) = source_size(i_source) / PI * 180
           obj(i_source) = comb_obj_init(comb_tmplmap_texture, nside, &
                amplitude(i_source), dilation=1.0, &
                alpha=alpha(i_source), beta=beta(i_source), gamma=gamma(i_source), &
                name=trim(comb_type), param=texture_params)
        end do

        
  end select

  ! Generate csky.
  if(beam_status) then
     if (comb_type == COMB_TYPE_OPT(7) .or. comb_type == COMB_TYPE_OPT(8)) then
        csky = comb_csky_init(obj, beam=beam)
     else
        csky = comb_csky_init(obj_mother, n_source, amplitude, &
             dilation_array, alpha, beta, gamma, beam=beam)
     end if
  else
     if (comb_type == COMB_TYPE_OPT(7) .or. comb_type == COMB_TYPE_OPT(8)) then
        csky = comb_csky_init(obj)
     else
        csky = comb_csky_init(obj_mother, n_source, amplitude, &
             dilation_array, alpha, beta, gamma)
     end if
  end if


  !---------------------------------------
  ! Generate output files
  !---------------------------------------

  call comb_csky_write_sky_obj(csky, filename_out)
  if (plot_all) then
     call comb_csky_write_sky_objs(csky, filename_out(1:len(trim(filename_out))-5))
  end if


  !---------------------------------------
  ! Free memory
  !---------------------------------------

  call comb_csky_free(csky)
  if (comb_type /= COMB_TYPE_OPT(7) .and. comb_type /= COMB_TYPE_OPT(8)) then
     call comb_obj_free(obj_mother)
  end if
  if(beam_status) call s2_pl_free(beam)
  deallocate(amplitude)
  deallocate(dilation_array)
  deallocate(alpha, beta, gamma)
  deallocate(source_size)
  if (allocated(obj)) deallocate(obj)


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
            write(*,'(a)') '                   [-include_size include_size (optional)]'
            write(*,'(a)') '                   [-plot_all plot_all (optional)]'
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

          case ('-include_size')
            read(arg,*) include_size

          case ('-plot_all')
            read(arg,*) plot_all

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program comb_objgen

