!------------------------------------------------------------------------------
!>  This is a template for a GKW diagnostic module.
!>
!> Before writing a new diagnostic or module, you should in
!> general try to find the most similar existing diagnostic in terms
!> of quantity calculated.  If at all possible,
!> you should generalize an existing diagnostic for your needs,
!> rather than copy/pasting and adapting it.
!>
!> To write a new diagnostic, first read some of the other
!> diagnostics, to learn how things are done.
!> Then fill the subroutines of this template with life and hook them
!> into the appropriate subroutines of diagnostic.F90 . You may delete
!> unneeded routines of this skeleton and not hook them into the
!> diagnostic.F90 . All routines are 'collective' in the MPI sense.
!> If a block has to be executed by only a single processor, you have
!> to write something like
!>    if(root_processor) then ... end if
!> or
!>    if(proc_subset(0,1,1,1,gsp(is))) then ... end if
!> Be aware of the parallelization in particular when you want to write
!> output data.
!>
!> Please try to encapsulate new diagnostics as much as possible, so
!> that they can be easily understood by others.
!> This means that you should not use public variables/routines of the
!> other diagnos_* modules, except from diagnos_generic.
!> Routines and switches which are used by more than one diagnos_* module
!> (e.g. generally useful output routines) should be placed
!> into diagnos_generic.
!>
!> It is really important put comments, *immediately*, at the time of
!> writing the code, to save other readers (and yourself) precious time
!> when they need to understand your code. Others should not have to waste
!> time thinking about ambiguities or non-trivial algorithms.
!>
!> The lifecycle of a module diagnos_foo during a run of the explicit
!> solver is usually
!>   1.  The namelist items provided by diagnos_foo are initialised
!>       with default values.
!>   2.  The diagnostic namelist is parsed and the default values of namelist
!>       items are overwritten.
!>   3.  bcast() is called.
!>   4.  check() is called.
!>   5.  allocate_mem() is called.
!>   6.  read_last_data() is called.
!>   7.  init() is called.
!>   8.  initial_output() is called.
!>   9.  output() is called repeatedly.
!>    evtl. screen_output() is called repeatedly.
!>   10. final_output() is called.
!>   11. finalize() is called.
!>
!> To make naming easier, the subroutine-names do not have a prefix
!> or suffix.  In order to call them from the diagnostic module und
!> avoid name clashes, you should use those subroutines like this:
!>
!> use diagnos_template, only : init_template => init
!> 
!> TESTCASES:
!>
!> Every diagnostic should be tested by one or several testcases.
!> In the GKW standard test set there are a few testcases which have
!> many diagnostics enabled. Please switch on your diagnostic and
!> add reference files when you think it is ready.
!>
!>------------------------------------------------------------------------------
!>
!> TEMPLATE MODULE DESCRIPTION
!> 
!> DESCRIPTION
!> 
!> [ Put a description of the diagnostic here. ]
!> [ You may also want to give a reference to the work it appears in. ]
!>
!> FURTHER NOTES: 
!>
!> [ In case your diagnostic really wants to use data calculated by
!>   other diagnostics (like fluxes), please mention it. ]
!>
!> LIMITATIONS: 
!>
!> [ Maybe you want to tell others what this module cannot (yet) do ]
!>
!------------------------------------------------------------------------------
module diagnos_template

  implicit none

  private
  
  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output, final_output, read_last_data
  public :: calc_smallstep
  public :: output, screen_output

  !> The general on/off switch for this diagnostic
  logical, save, public :: lcalc_template

  !> The logical unit number, used to output data
  integer, save :: lun_template

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    ! Usually the diagnostic should be switched
    ! off by default.

    lcalc_template = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_template,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_warn

    if(.not. lcalc_template) return

    if(.true. .eqv. .false.) then
      call gkw_warn('The template diagnostic does not work with xyz, switched &
         & off.')
      ! switch yourself off
      lcalc_template = .false.
    end if

    ! Typical things to warn about:
    ! Does it work with nonlin. runs?
    ! Does it work with nonspectral/spectral runs?
    ! Does it work with normalized runs?
    ! Does it work with restarted runs?
    ! Does it work with modebox=F runs?
    ! Does it work with electromagnetic runs?
    ! Does it work with trace species (where de(:,:) = 0)?

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid, lwrite_output1
    use global, only : DISTRIBUTION
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, X_GHOSTCELLS
    logical, intent(inout) :: requirements(:,:)

    if(.not. lcalc_template) return

    call register_tag_range(the_number_of_tags_this_diagnostic_needs, &
       & tag_range_start, tag_range_end_inkl)

    ! specify what your diagnostic needs. This is important in
    ! particular to provide it with ghost cell data if it computes
    ! derivatives. Important: Only set things to .true., do not set
    ! them to .false. if they are not needed by this diagnostic.
    ! Example:
    requirements(DISTRIBUTION,LOCAL_DATA) = .true.
    requirements(PHI_FIELD,LOCAL_DATA) = .true.
    ! requirements(APAR_FIELD,LOCAL_DATA) = .true.
    ! requirements(BPAR_FIELD,LOCAL_DATA) = .true.
    ! requirements(F_DISTRIBUTION,S_GHOSTCELLS) = .true.
    requirements(PHI_FIELD,S_GHOSTCELLS) = .true.
    ! requirements(APAR_FIELD,S_GHOSTCELLS) = .true.
    ! requirements(BPAR_FIELD,S_GHOSTCELLS) = .true.
    ! requirements(F_DISTRIBUTION,X_GHOSTCELLS) = .true.
    !...
    
    if(root_processor) then
      ! Open logical units here.

      ! For per-timestep quantities, open the logical unit only if
      ! lwrite_output1 is set.
      
      ! A PROPOSAL FOR NAMING RULES OF NEW LOGICAL UNITS: Note that the logical
      ! unit name does not
      ! contain the .dat, .bin or any other suffix.  Moreover the
      ! logical unit name should not contain dots (because in
      ! scripting languages the data is likely to be stored in a
      ! structure and there the dot is an operator and needs to be converted).
      ! Good rule:
      ! choose a logical unit name which would be a valid Fortran
      ! variable name, and if in future there might also be other logical
      ! units of the same quantity, you should also specify the
      ! representation (a combination of x,y,s,sp,kx,ky,ks,vp,m) in
      ! a short notation, if it is possible. An enumerating counter should be
      ! inserted right after the shortname.
      ! 
      !   <shortname><filecount>_<repr>
      !  or <shortname>_sp[0-9]*_<repr>  (with a species index)
      !  or <shortname>_s[0-9]*_<repr>   (with a coordinate index)
      !  ...
      !
      ! e.g. entr_ky, entr_x, phisq_kykx, dens_xy, temp_ykx,
      !      foo0001, bar001234_xy, dens_sp01_s24_xy, distr009999_sp01_vpm
      
      call open_real_lu('template', 'diagnostic/diagnos_template', &
         & (/ nmod, n_x_grid /), &
         & ascii_fmt, lun_template)
      ! if its a 1D quantity:
      ! (replace with the right grid names, of course)
      !call attach_metadata_grid(lun_template, 'time', ascii_fmt)
      ! if its a 2D quantity:
      call attach_metadata_grid(lun_template, 'time', 'krho', ascii_fmt)
      ! if its a 3D quantity:
      !call attach_metadata_grid(lun_template, 'time', 'krho', 'kxrh', ascii_fmt)
      
      ! the unit in terms of reference quantities, should be valid LaTeX:
      call attach_metadata(lun_template, phys_unit_key, &
         & 'v_{th,ref}n_{R,0}R_{ref}B_{ref}', ascii_fmt)
      call attach_metadata(lun_template, description_key, &
         & 'Not tremendously long description of the output quantity here.', &
         & ascii_fmt)
      ! and some (more informal) comments there, if you want:
      call attach_metadata(lun_template, comments_key, &
         & not_avail, ascii_fmt)
      ! and you may attach even more infos:
      call attach_metadata(lun_template, 'important_parameter', &
         & 123.456, ascii_fmt)
    end if
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    integer :: ierr

    if(.not. lcalc_template) return
    
    ! allocate(local_array_name(nmod,number_of_species),stat=ierr)
    ! if (ierr /= 0) &
    !    & call gkw_abort('diagnos_template :: could not allocate array')
  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine read_last_data()

    if(.not. lcalc_template) return
    
  end subroutine read_last_data

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, ascii_fmt

    if(.not. lcalc_template) return
    
    ! deallocate all arrays of this diagnostic
    !if(allocated(local_array_name)) deallocate(local_array_name)

    ! be nice and close all logical units
    call close_lu(lun_template, ascii_fmt)
  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

    if(.not. lcalc_template) return

  end subroutine initial_output

  !--------------------------------------------------------------------
  !> In contrast to the subroutine output(), the argument of
  !> the final_output() subroutine is optional.
  !> The number is provided if the code wants to print the result of
  !> an eigenmode calculation.
  !> If the GKW run is explicit (or implicit) then the number argument
  !> is not provided.
  !--------------------------------------------------------------------
  subroutine final_output(number)

    integer, intent(in), optional :: number

    if(.not. lcalc_template) return

    if (present(number)) then
      ! output for every eigenmode
    else
      ! one single output
    end if

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(i_smallstep)
    integer, intent(in) :: i_smallstep

    if(.not. lcalc_template) return
    
    ! if (i_smallstep == naverage - 1) then

    ! else if (i_smallstep == naverage) then

    ! end if

  end subroutine calc_smallstep

  !--------------------------------------------------------------------
  !> The routine calc_largestep() is to be called repeatedly, after every
  !> large timestep, and here the diagnostic should calculate its
  !> quantities.
  !> 
  !> Splitting calculation and output is to be encouraged, particularly when
  !>  - the output is higher dimensional, and requires MPI or MPI-IO
  !>  - same data is output in different dimensional slices
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use dist, only : phi, apar, bpar
    
    if(.not. lcalc_template) return

    ! phi, apar, bpar are already filled, just read them, do not change them

  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()

    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpiinterface, only : root_processor

    if(.not. lcalc_template) return

    if(root_processor) then
      call append_chunk(lun_template, local_array_name, xy_fmt, ascii_fmt)
    end if

  end subroutine output

  !--------------------------------------------------------------------
  !> If useful, print something to the stdout.
  !--------------------------------------------------------------------
  subroutine screen_output()
    use mpiinterface, only : root_processor

    if(.not. lcalc_template) return
    
    if( root_processor) then
      !write(*,*)'Quantity X:', X
    end if
  end subroutine screen_output


end module diagnos_template
