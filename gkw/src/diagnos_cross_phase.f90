!------------------------------------------------------------------------------
!> Calculate probability density functions (PDFs) for cross-phases
!>
!> Consider two fluctuations
!>
!>   A = exp[i alpha]
!>   B = exp[i beta]
!>
!> Their phase-difference alpha-beta ("cross-phase") can be computed as
!>
!>   alpha - beta = tan^{-1}[ Im{exp[i(alpha-beta)]} / Re{exp[i(alpha-beta)]} ]
!>                = tan^{-1}[ Im{A/B}/Re{A/B} ]
!>
!> A weighted probability density function of the cross-phase is computed, where
!> the weight of each contribution is determined from the strength of the
!> two fluctuating quantities involved.
!>
!> Cross-phases are sometimes measured in the literature to assess how "linear"
!> certain modes are.
!>
!> Cross-phases between electron density and temperature can be measured
!> experimentally. Therefore this diagnostic can be useful for validation
!> purposes.
!>
!> The PDF is continuously updated at every large timestep and finally
!> output at the end of the run.
!>
!> To compute momenta consistently, this diagnostic is connected to
!> diagnos_moments.
!> For weighting contributions to the PDF in linear runs, this diagnostic may
!> use growth rates computed previously by diagnos_growth_freq .
!>
!> Random literature examples which use this kind of measurement are, e.g.
!>   T. Dannert, F. Jenko, Phys. Plasmas 12, 072309 (2005)
!>   A. Navarro et al., Phys. Plasmas 22, 042513 (2015)
!>
!------------------------------------------------------------------------------
module diagnos_cross_phase
  use diagnos_generic, only : max_token_length

  implicit none

  private
  
  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output, final_output, read_last_data
  public :: output
  public :: get_field_or_moment

  integer, parameter :: max_cross_phases = 32

  ! this hardcoded switch can be changed to deactivate weighting
  logical, parameter :: weigh_pdf_across_y = .true.
  ! this hardcoded switch can be changed to deactivate weighting
  logical, parameter :: weigh_pdf_along_s_x = .true.
  integer, parameter :: WEIGHTED_SX = 2
  integer, parameter :: WEIGHTED_YSX = 1

  ! naturally the zonal flow sucks up free energy and thus most of the
  ! weight if weigh_pdf_across_y=true . Use this factor here to ignore
  ! the ZF (i.e. do not measure it's cross-phase) by setting this to
  ! 0.0 . Do not do anything special to the ZF by setting this to 1.0
  ! .
  real, parameter :: special_zf_weight_factor = 1.0

  
  !> A namelist parameter: how many classes for the phase should the PDF have
  integer, save, public :: cross_phase_nclasses

  ! A namelist parameter: output timetrace of cross phases
  logical, save, public :: cross_phase_timetrace

  !> A namelist parameter: list, which indicates the cross-phases that are to be
  !> calculated. This makes use of the indices PHI_FIELD, etc.  It
  !> would be more readable in the input file, to give strings, but it
  !> is way easier to implement a list of integers.
  character(len=max_token_length), save, &
     & public :: cross_phases(2*max_cross_phases)
  integer, save, public :: cross_phases_parsed(2,max_cross_phases)

  !> this is the effective number of cross phases to compute
  !> (respecting that there are number_of_species densities,
  !> temperatures, etc. )
  integer, save :: n_cross_phases

  !> the buffer to compute the pdf of the phase value for each nmod
  !> and each field/moment-combination
  real, save, allocatable :: pdf(:,:,:,:)

  !> logical unit numbers to save the PDFs of each timestep
  integer, save, allocatable :: timetrace_luns(:)

  !> this is global in species, but local in space dimensions
  complex, save, allocatable :: A_buf(:, :, :, :)
  !> this is global in species, but local in space dimensions
  complex, save, allocatable :: B_buf(:, :, :, :)

  !> A datatype that represents a (nsp,nmod,ns,nx) block in a
  !> (number_of_species,nmod,ns,nx) array.
  integer, save :: mpi_dtype

contains

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

    ! is switched off by default
    cross_phases = ''
    cross_phase_nclasses = 181

    cross_phase_timetrace = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    integer :: i
    do i = 1, 2*max_cross_phases
      call mpibcast(cross_phases(i),max_token_length)
    end do
    call mpibcast(cross_phase_nclasses, 1)
    call mpibcast(cross_phase_timetrace, 1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use global, only : root_and_verbose
    use general, only : gkw_warn, gkw_abort
    use diagnos_generic, only : START_MOMENTS
    
    use diagnos_generic, only : DENSITY_MOMENT, TEMPERATURE_MOMENT
    use diagnos_generic, only : PAR_TEMPERATURE_MOMENT, PERP_TEMPERATURE_MOMENT
    use diagnos_generic, only : CURRENT_MOMENT
    use diagnos_generic, only : PASSING_DENSITY_MOMENT, TRAPPED_DENSITY_MOMENT
    use diagnos_generic, only : DENSITY_GA_MOMENT, DENSITY_POLAR_MOMENT
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    
    use diagnos_generic, only : DENSITY_MOMENT_TOKEN, TEMPERATURE_MOMENT_TOKEN
    use diagnos_generic, only : PAR_TEMPERATURE_MOMENT_TOKEN, PERP_TEMPERATURE_MOMENT_TOKEN
    use diagnos_generic, only : CURRENT_MOMENT_TOKEN
    use diagnos_generic, only : PASSING_DENSITY_MOMENT_TOKEN, TRAPPED_DENSITY_MOMENT_TOKEN
    use diagnos_generic, only : DENSITY_GA_MOMENT_TOKEN, DENSITY_POLAR_MOMENT_TOKEN
    use diagnos_generic, only : PHI_FIELD_TOKEN, APAR_FIELD_TOKEN, BPAR_FIELD_TOKEN
    use diagnos_generic, only : PHI_GA_FIELD_TOKEN, APAR_GA_FIELD_TOKEN, BPAR_GA_FIELD_TOKEN
    use diagnos_growth_freq, only : lgrowth_rates
    use control, only : nlapar, nlbpar, spectral_radius, non_linear
    use restart, only : restarted
    use grid, only : number_of_species, nmod
    integer :: i, n, m

    n_cross_phases = 0
    do i = 1, max_cross_phases
      m = 1
      do n = 1,2
        select case(cross_phases(2*(i-1) + n))
        case(PHI_FIELD_TOKEN)
          cross_phases_parsed(n,i) = PHI_FIELD
          m = m * 1
        case(APAR_FIELD_TOKEN)
          if(nlapar) then
            cross_phases_parsed(n,i) = APAR_FIELD
            m = m * 1
          else
            call gkw_abort('cross_phases: run does not compute '//APAR_FIELD_TOKEN)
          end if
        case(BPAR_FIELD_TOKEN)
          if(nlbpar) then
            cross_phases_parsed(n,i) = BPAR_FIELD
            m = m * 1
          else
            call gkw_abort('cross_phases: run does not compute '//BPAR_FIELD_TOKEN)
          end if
        case(PHI_GA_FIELD_TOKEN)
          call gkw_abort('cross_phases: not implemented and probably does not &
             & make much sense')
        case(APAR_GA_FIELD_TOKEN)
          call gkw_abort('cross_phases: not implemented and probably does not &
             & make much sense')
        case(BPAR_GA_FIELD_TOKEN)
          call gkw_abort('cross_phases: not implemented and probably does not &
             & make much sense')
        case('')
          cross_phases_parsed(n,i) = 0
        case(DENSITY_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + DENSITY_MOMENT
          m = m * number_of_species
        case(TEMPERATURE_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + TEMPERATURE_MOMENT
          m = m * number_of_species
        case(PAR_TEMPERATURE_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + PAR_TEMPERATURE_MOMENT
          m = m * number_of_species
        case(PERP_TEMPERATURE_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + PERP_TEMPERATURE_MOMENT
          m = m * number_of_species
        case(CURRENT_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + CURRENT_MOMENT
          m = m * number_of_species
        case(PASSING_DENSITY_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + PASSING_DENSITY_MOMENT
          m = m * number_of_species
        case(TRAPPED_DENSITY_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + TRAPPED_DENSITY_MOMENT
          m = m * number_of_species
        case(DENSITY_GA_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + DENSITY_GA_MOMENT
          m = m * number_of_species
          if(.not. spectral_radius) then
            call gkw_abort('cross_phases: dens_ga is not implemented for nonspectral')
          end if
        case(DENSITY_POLAR_MOMENT_TOKEN)
          cross_phases_parsed(n,i) = START_MOMENTS + DENSITY_POLAR_MOMENT
          m = m * number_of_species
          if(.not. spectral_radius) then
            call gkw_abort('cross_phases: dens_polar is not implemented for nonspectral')
          end if
        case default
          call gkw_abort('cross_phases: wrong token given')
        end select
      end do
      if(cross_phases_parsed(1,i) == 0 .neqv. cross_phases_parsed(2,i) == 0) then
        call gkw_abort('cross_phases: cannot group tokens in pairs')
      end if
      if(all(cross_phases_parsed(:,i) /= 0)) then
        n_cross_phases = n_cross_phases + m
        if(root_and_verbose) then
          write(*,*) 'Cross-phase diagnostic: will compute cross-phase of ', &
             & cross_phases(2*(i-1)+1), cross_phases(2*(i-1)+2)
        end if
      end if
    end do

    if(n_cross_phases > 0 .and. restarted) then
      call gkw_warn('cross_phases: This diagnostics computes a PDF which is &
         & output only at the end of the run. As this run is restarted, you &
         & may end up with more than 1 output dataset.')
    end if

    if(n_cross_phases > 0 .and. mod(cross_phase_nclasses,2) /= 1) then
      call gkw_abort('cross_phases: cross_phase_nclasses must be odd.')
    end if

    if(weigh_pdf_across_y .and. .not.non_linear .and. nmod > 1) then
      call gkw_warn('cross_phases: This diagnostic enforces lgrowth_rates=T &
         & for weighting the pdf')
      lgrowth_rates = .true.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : MPICOMPLEX_X
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use diagnos_generic, only : LOCAL_DATA
    use global, only : id_sp, id_mod, id_x, id_s
    use mpidatatypes, only : create_subarray_datatype
    use control, only : non_linear
    use grid, only : nx, ns, proc_subset
    
    logical, intent(inout) :: requirements(:,:)

    if(n_cross_phases == 0) return

    ! specify what your diagnostic needs. This is important in
    ! particular to provide it with ghost cell data if it computes
    ! derivatives. Important: Only set things to .true., do not set
    ! them to .false. if they are not needed by this diagnostic.
    ! Example:
    if(any(cross_phases_parsed == PHI_FIELD)) then
      requirements(PHI_FIELD,LOCAL_DATA) = .true.
    end if
    if(any(cross_phases_parsed == APAR_FIELD)) then
      requirements(APAR_FIELD,LOCAL_DATA) = .true.
    end if
    if(any(cross_phases_parsed == BPAR_FIELD)) then
      requirements(BPAR_FIELD,LOCAL_DATA) = .true.
    end if

    if(.true. .or. non_linear) pdf = 0.0

    if(proc_subset(0,0,1,1,0)) then
      ! create subarray MPI datatype for 4d spectral field
      call create_subarray_datatype(MPICOMPLEX_X, mpi_dtype, &
         & id_sp, id_mod, id_s, id_x, &
         & global_and_local_xsize=nx, global_and_local_ssize=ns)
    end if

    if(cross_phase_timetrace) then
      ! it is easier to just initialise these luns with a zero (to
      ! mark them as not yet opened) and actually open them in output_pdf()
      timetrace_luns = 0
    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general, only : gkw_abort
    use grid, only : nmod, ns, nx, number_of_species, proc_subset
    integer :: ierr

    if(n_cross_phases == 0) return

    ! the 0 index part is a temporary buffer
    allocate(pdf(2,cross_phase_nclasses,nmod,0:n_cross_phases),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_cross_phase :: could not allocate pdf')

    if(cross_phase_timetrace) then
      allocate(timetrace_luns(n_cross_phases),stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_cross_phase :: could not allocate &
         &timetrace_luns')
    end if

    if(proc_subset(0,0,1,1,1)) then
      allocate(A_buf(number_of_species, nmod, ns, nx),stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_cross_phase :: could not allocate A_buf')
      allocate(B_buf(number_of_species, nmod, ns, nx),stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_cross_phase :: could not allocate B_buf')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine read_last_data()

    if(n_cross_phases == 0) return

    ! FIXME this might evtl even make sense to implement for this
    ! diagnostic: otherwise one cannot continue a run with this
    ! diagnostic and get the same data as from a single run.
    
  end subroutine read_last_data

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, ascii_fmt
    integer :: m

    if(n_cross_phases == 0) return
    
    ! deallocate all arrays of this diagnostic
    if(allocated(pdf)) deallocate(pdf)

    if(cross_phase_timetrace) then
      do m = 1, size(timetrace_luns)
        call close_lu(timetrace_luns(m),ascii_fmt)
      end do
      if(allocated(timetrace_luns)) deallocate(timetrace_luns)
    end if
    
    if(allocated(A_buf)) deallocate(A_buf)
    if(allocated(B_buf)) deallocate(B_buf)

  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

    if(n_cross_phases == 0) return

    call output_degr_grid()
    
  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output(number)
    use control, only : non_linear
    integer, intent(in), optional :: number

    if(n_cross_phases == 0) return
    
    if(.not. (.true. .or. non_linear)) then
      pdf = 0
      call calc_largestep(0.0)
      if (present(number)) then
        ! output for every eigenmode
        call output_all_pdfs(number)
      else
        ! one single output
        call output_all_pdfs()
      end if
    else
      ! calculation of the pdf is done at every timestep
      call output_all_pdfs()
    end if

  end subroutine final_output

  !--------------------------------------------------------------------
  !> 
  !--------------------------------------------------------------------
  subroutine calc_largestep(last_time)
    use global, only : r_tiny
    use grid, only : number_of_species, nmod, nx, ns, proc_subset
    use grid, only : n_x_grid, n_s_grid
    use control, only : time, non_linear
    use global, only : BPAR_GA_FIELD
    use diagnos_generic, only : tokens
    use constants, only : pi
    use global, only : int2char_zeros, int2char, r_tiny, root_and_verbose
    use general, only : gkw_warn
    use mpiinterface, only : mpiallreduce_sum, mpireduce_sum_inplace
    use mpicomms, only : COMM_S_NE_X_NE
    use diagnos_growth_freq, only : growth_rates
    use mode, only : iyzero, ixzero, mode_label_G
    !> should be either zero or control::last_largestep_time
    real, intent(in) :: last_time
    integer :: i, n, class
    integer :: m, m1, m2, im1, im2
    integer :: imod, ix, is
    real :: phase
    real :: dprobability, s_x_weight = 0.0, s_x_weight_local_total(nmod)
    real :: s_x_weight_norm(nmod), y_weight(nmod)
    real :: delta_angle
    integer :: fpe_counter
    fpe_counter = 0
    
    if(n_cross_phases == 0) return

    ! phi, apar, bpar are already filled, just read them, do not change them
    n = 0
    all_token_pairs: do i = 1, max_cross_phases
      valid_pair: if(cross_phases_parsed(1,i) /= 0) then
        if(cross_phases_parsed(1,i) <= BPAR_GA_FIELD) then
           m1 = 1
        else
          m1 = number_of_species
        end if
        if(cross_phases_parsed(2,i) <= BPAR_GA_FIELD) then
           m2 = 1
        else
          m2 = number_of_species
        end if
        species_loop1: do im1 = 1, m1
          ! put the moment or field into the A buffer
          !A_buf = get_field_or_moment(cross_phases_parsed(1,i))
          call get_field_or_moment(cross_phases_parsed(1,i), A_buf)
          
          species_loop2: do im2 = 1, m2
            m = n+im1+(im2-1)*m2

            if(root_and_verbose) then
              ! debug info
              write(*,*) 'Computing cross-phase #'//int2char_zeros(m,2)//' of ', &
                 & tokens(cross_phases_parsed(1,i)), 'species', im1, &
                 & tokens(cross_phases_parsed(2,i)), 'species', im2
            end if

            ! put the moment or field into the B buffer
            !B_buf = get_field_or_moment(cross_phases_parsed(2,i))
            call get_field_or_moment(cross_phases_parsed(2,i), B_buf)

            ! while all processes did go into the computation of the
            ! moments, now they should be available on the procs
            ! working on global species 1.

            if(proc_subset(0,0,1,1,1)) then
              if(weigh_pdf_across_y) then
                ! weighted probability: the weight factor is the
                ! product of the strength of the two considered
                ! fields. This is useful for nonlinear runs.
                if(non_linear) then
                  do imod = 1, nmod
                    call mpiallreduce_sum(sum(abs(B_buf(im2,imod,:,:)* &
                       & A_buf(im1,imod,:,:))),y_weight(imod), 1, COMM_S_NE_X_NE)
                  end do
                  ! allow to ignore the ZF, as it may become super strong.
                  if(iyzero /= 0) then
                    y_weight(iyzero) = y_weight(iyzero) * special_zf_weight_factor
                  end if
                else if(nmod > 1) then
                  ! linear runs are typically normalized. Instead of
                  ! with the strength of the fluctuations, it may be
                  ! good to weigh using the growth rate.
                  do imod = 1, nmod
                    y_weight(imod) = max(growth_rates(mode_label_G(imod,ixzero)),&
                       & 0.0)
                  end do
                else
                  y_weight = 1.0
                end if

                if(abs(sum(y_weight)) > r_tiny) then
                  y_weight = y_weight / sum(y_weight)
                else
                  y_weight = 1.0/nmod
                end if
              else
                ! equal probability for all toroidal modes
                y_weight = 1.0/nmod
              end if

              if(weigh_pdf_along_s_x) then
                do imod = 1, nmod
                  s_x_weight_local_total(imod) = sum(abs(B_buf(im2,imod,:,:)* &
                     & A_buf(im1,imod,:,:)))
                end do
                call mpiallreduce_sum(s_x_weight_local_total,s_x_weight_norm, &
                   & nmod, COMM_S_NE_X_NE)
              end if

              ! calculate the probability density for the current timestep,
              ! store it into pdf(:,:,:,0)
              pdf(:,:,:,0) = 0.0
              sloop: do is = 1, ns
                xloop: do ix = 1, nx
                  yloop: do imod = 1, nmod
                    if(weigh_pdf_along_s_x) then
                      ! weighted probability: the weight factor is the
                      ! product of the two considered fields.
                      if(abs(s_x_weight_norm(imod)) > r_tiny) then
                        ! this divisor can be zero, e.g. in a linear
                        ! simulation with nmod>1
                        s_x_weight = abs(B_buf(im2,imod,is,ix) * &
                           & A_buf(im1,imod,is,ix)) / s_x_weight_norm(imod)
                      end if
                    else
                      ! equal probability at all s,x points for this imod
                      s_x_weight = 1.0/(n_x_grid*n_s_grid)
                    end if

                    if(abs(B_buf(im2,imod,is, ix)) > r_tiny) then
                      if(abs(real(A_buf(im1,imod,is,ix)/&
                         & B_buf(im2,imod,is,ix))) > r_tiny) then

                        phase = atan2(aimag(&
                           & A_buf(im1,imod,is,ix)/B_buf(im2,imod,is,ix)), &
                           & real(A_buf(im1,imod,is,ix)/B_buf(im2,imod,is,ix)))
                        ! the atan returns values between -pi/2 and pi/2
                        ! but the atan2 between -pi and pi

                        ! define classes to make a histogram.
                        delta_angle = 2*pi / (cross_phase_nclasses-1)

                        ! the classes go from 1 to cross_phase_nclasses
                        class = ceiling((phase + delta_angle*0.5)/delta_angle) + &
                           & int((cross_phase_nclasses-1)*0.5)

                        ! the outmost classes have just half the width.
                        if(class == 1 .or. class == cross_phase_nclasses) then
                          delta_angle = delta_angle * 0.5
                        end if

                        dprobability = s_x_weight*y_weight(imod)
                        pdf(WEIGHTED_YSX,class,imod,0) = pdf(WEIGHTED_YSX,class,imod,0) + &
                           & dprobability/delta_angle
                        
                        dprobability = s_x_weight
                        pdf(WEIGHTED_SX,class,imod,0) = pdf(WEIGHTED_SX,class,imod,0) + &
                           & dprobability/delta_angle

                      else
                        fpe_counter = fpe_counter + 1
                        if(root_and_verbose) then
                          write(*,*) "Cross-phase diagnostic: fpe1 at", &
                             & is, ix, imod
                        end if
                        ! the expression
                        ! for the phase cannot be computed because of
                        ! division by zero
                      end if
                    else
                      fpe_counter = fpe_counter + 1
                      if(root_and_verbose) then
                        write(*,*) "Cross-phase diagnostic: fpe2 at", &
                           & is, ix, imod
                      end if
                    end if
                  end do yloop
                end do xloop
              end do sloop
              
              if(cross_phase_timetrace .and. proc_subset(0,0,1,1,1)) then
                ! combine the local PDFs from different parts of the
                ! grid by summing
                call mpireduce_sum_inplace(pdf(:,:,:,0), &
                   & shape(pdf(:,:,:,0)), COMM_S_NE_X_NE)

                ! output the current PDF
                call output_pdf(pdf(WEIGHTED_YSX,:,:,0),i,m,m1,im1,m2,im2, .true.,WEIGHTED_YSX)
                call output_pdf(pdf(WEIGHTED_SX,:,:,0),i,m,m1,im1,m2,im2, .true.,WEIGHTED_SX)
              end if

              ! combine the current PDF with the running PDF
              pdf(:,:,:,m) = (pdf(:,:,:,m) * last_time + &
                 & pdf(:,:,:,0) * (time-last_time)) / time

            end if
          end do species_loop2
        end do species_loop1
        ! the pair-counter to access the right part of pdf(:,:,:,:)
        n = n + m1*m2

      end if valid_pair
    end do all_token_pairs

    if(fpe_counter > 0) then
      call gkw_warn('Cross-phase diagnostic: division-by-zero count: ' &
         & //int2char(fpe_counter))
    end if

  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !> For debugging: check if the PDF integrates to 1
  !--------------------------------------------------------------------
  subroutine check_integrates_to_one(weight_mode)
    use global, only : int2char_zeros
    use constants, only : pi
    use grid, only : nmod
    integer, intent(in) :: weight_mode
    integer :: m, i
    real :: s(nmod)
    do m = 1, ubound(pdf,4)
      s = 0
        do i = 1, cross_phase_nclasses
          if(i == 1 .or. i == cross_phase_nclasses) then
            s = s + pdf(weight_mode,i,:,m) * (0.5 * 2*pi / (cross_phase_nclasses-1))
          else
            s = s + pdf(weight_mode,i,:,m) * (2*pi / (cross_phase_nclasses-1) )
          end if
      end do
      write(*,*) "Cross-phase diagnostic: Weight-mode "//&
         & int2char_zeros(weight_mode,1)&
         & //", sum of pdf #"//int2char_zeros(m,2), &
         & " is ", sum(s), "or", s
    end do
  end subroutine check_integrates_to_one

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine get_field_or_moment(iquantity,ret_buf)
    use general, only : gkw_abort
    use diagnos_generic, only : START_MOMENTS
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use diagnos_moments, only : get_4d_moment
    use mpiinterface, only : mpireduce_sum_inplace, gather_array
    use mpicomms, only : COMM_VPAR_NE_MU_NE, COMM_SP_NE
    use grid, only : nsp, nmod, ns, nx, proc_subset
    use dist, only : phi, apar, bpar
    integer, intent(in) :: iquantity
    complex :: local_array(nsp,nmod,ns,nx)
    !> ret_buf is global in the species, but local in the space dimensions.
    !> this is only returned on proc_subset(0,0,1,1,1)
    !complex :: ret_buf(number_of_species, nmod, ns, nx)
    complex, intent(out) :: ret_buf(:,:,:,:)
    integer :: i, ix, imod

    select case(iquantity)
    case(PHI_FIELD)
      ! I cannot simply take phi(nmod,nx,ns) from dist because the moments
      ! are computed on (nmod,ns,nx) arrays, thanks to history.
      if(proc_subset(0,0,1,1,1)) then
        do i = 1, ns
          do ix = 1, nx
            do imod = 1, nmod
              ret_buf(1,imod,i,ix) = phi(imod,ix,i)
            end do
          end do
        end do
      end if
    case(APAR_FIELD)
      if(proc_subset(0,0,1,1,1)) then
        do i = 1, ns
          do ix = 1, nx
            do imod = 1, nmod
              ret_buf(1,imod,i,ix) = apar(imod,ix,i)
            end do
          end do
        end do
      end if
    case(BPAR_FIELD)
      if(proc_subset(0,0,1,1,1)) then
        do i = 1, ns
          do ix = 1, nx
            do imod = 1, nmod
              ret_buf(1,imod,i,ix) = bpar(imod,ix,i)
            end do
          end do
        end do
      end if
    case(PHI_GA_FIELD)
      call gkw_abort('cross_phases: not implemented and probably does not make&
         & much sense')
    case(APAR_GA_FIELD)
      call gkw_abort('cross_phases: not implemented and probably does not make&
         & much sense')
    case(BPAR_GA_FIELD)
      call gkw_abort('cross_phases: not implemented and probably does not make&
         & much sense')
    case default
      local_array = get_4d_moment(iquantity - START_MOMENTS)
      call mpireduce_sum_inplace(local_array, shape(local_array), &
         & COMM_VPAR_NE_MU_NE)
      if(proc_subset(0,0,1,1,0)) then
        call gather_array(ret_buf, local_array, mpi_dtype, COMM_SP_NE, &
           & .false., .true.)
        ! now the procs proc_subset(0,0,1,1,1) should have an array
        ! which is global in the species, but local in space.
      end if
      
    end select

  end subroutine get_field_or_moment

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()
    use control, only : non_linear, last_largestep_time

    if(n_cross_phases == 0) return

    if(.true. .or. non_linear) then
      call calc_largestep(last_largestep_time)
    end if

  end subroutine output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_all_pdfs(number)
    use mpiinterface, only : mpireduce_sum_inplace
    use mpicomms, only : COMM_S_NE_X_NE
    use global, only : BPAR_GA_FIELD
    use grid, only : number_of_species, proc_subset
    use global, only : int2char_zeros, root_and_verbose
    integer, intent(in), optional :: number

    integer :: i, n
    integer :: m, m1, m2, im1, im2

    if(.not.cross_phase_timetrace .and. proc_subset(0,0,1,1,1)) then
      ! combine the local PDFs from different parts of the
      ! grid by summing
      call mpireduce_sum_inplace(pdf(:,:,:,1:n_cross_phases), &
         & shape(pdf(:,:,:,1:n_cross_phases)), COMM_S_NE_X_NE)
    end if

    n = 0

    all_token_pairs: do i = 1, max_cross_phases
      valid_pair: if(cross_phases_parsed(1,i) /= 0) then
        if(cross_phases_parsed(1,i) <= BPAR_GA_FIELD) then
          m1 = 1
        else
          m1 = number_of_species
        end if
        if(cross_phases_parsed(2,i) <= BPAR_GA_FIELD) then
          m2 = 1
        else
          m2 = number_of_species
        end if
        species_loop1: do im1 = 1, m1

          species_loop2: do im2 = 1, m2

            m = n+im1+(im2-1)*m2

            if(present(number)) then
              call output_pdf(pdf(WEIGHTED_YSX,:,:,m),i,m,m1,im1,m2,im2,.false.,WEIGHTED_YSX,number)
              call output_pdf(pdf(WEIGHTED_SX,:,:,m),i,m,m1,im1,m2,im2,.false.,WEIGHTED_SX,number)
            else
              call output_pdf(pdf(WEIGHTED_YSX,:,:,m),i,m,m1,im1,m2,im2,.false.,WEIGHTED_YSX)
              call output_pdf(pdf(WEIGHTED_SX,:,:,m),i,m,m1,im1,m2,im2,.false.,WEIGHTED_SX)
            end if
          end do species_loop2
        end do species_loop1
        ! the pair-counter to access the right part of pdf(:,:,:,:)
        n = n + m1*m2

      end if valid_pair
    end do all_token_pairs

    if(root_and_verbose) then
      call check_integrates_to_one(WEIGHTED_YSX)
      call check_integrates_to_one(WEIGHTED_SX)
    end if

  end subroutine output_all_pdfs

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_pdf(data, i, m, m1, im1, m2, im2, do_append, weight_mode, &
     & number)
    use mpiinterface, only : root_processor
    use io, only : output_array, attach_metadata, open_real_lu, append_chunk
    use io, only : xy_fmt, ascii_fmt, description_key, comments_key, phys_unit_key
    use io, only : lu_exists
    use diagnos_generic, only : attach_metadata_grid, tokens
    use global, only : int2char_zeros
    use general, only : gkw_abort
    !> the pdf data
    real, intent(in) :: data(:,:)
    !> index of the cross_phase pair, as parsed
    integer, intent(in) :: i
    !> index of the pdf(:,:,:,m) that is to be output
    integer, intent(in) :: m
    !> m1 is 1 for fields and m1 > 1 if quantity exists for several species
    integer, intent(in) :: m1
    !> im1 is 1 for fields, and im1 >= 1 to indicate which species
    integer, intent(in) :: im1
    !> the same for the second quantity involved in the cross phase
    integer, intent(in) :: m2, im2
    !> append to or replace existing dataset
    logical, intent(in) :: do_append
    !> a parameter indicating which kind of weighting was applied
    integer, intent(in) :: weight_mode
    !> eigenmode number
    integer, intent(in), optional :: number
    
    character(len=11) :: sp1_snippet, sp2_snippet
    character(len=12) :: eim_number
    character(len=10) :: timetrace_snippet
    character(len=20) :: weight_snippet
    character(len=180) :: luname

    logical :: exists

    if(present(number)) then
      eim_number = '_eim'//int2char_zeros(number,1)
    else
      eim_number = ''
    end if

    if(m1 == 1) then
      sp1_snippet = ''
    else
      sp1_snippet = '_sp'//int2char_zeros(im1,2)
    end if
    
    if(m2 == 1) then
      sp2_snippet = ''
    else
      sp2_snippet = '_sp'//int2char_zeros(im2,2)
    end if

    if(do_append) then
      timetrace_snippet = 'timetrace_'
    else
      timetrace_snippet = ''
    end if

    if(weight_mode == WEIGHTED_YSX) then
      weight_snippet = '_ysx_weighted'
    else if(weight_mode == WEIGHTED_SX) then
      weight_snippet = '_sx_weighted'
    else
      call gkw_abort('cross_phases: unknown weight mode')
    end if

    if(root_processor) then
      luname = 'cross_phase_pdf_'//trim(timetrace_snippet)//&
         & trim(tokens(cross_phases_parsed(1,i)))//&
         & trim(sp1_snippet)//&
         & '_'// &
         & trim(tokens(cross_phases_parsed(2,i)))//&
         & trim(sp2_snippet)//&
         & trim(weight_snippet)//&
         & trim(eim_number)
      if(do_append) then
        if(.not. cross_phase_timetrace) then
          call gkw_abort('cross_phases: do_append but not cross_phase_timetrace?&
             & This was not intended.')
        end if
        if(timetrace_luns(m) == 0) then
          exists = .false.
          call open_real_lu(trim(luname),&
             & 'diagnostic/diagnos_cross_phase', &
             & shape(data), ascii_fmt, timetrace_luns(m))
        else
          exists = .true.
        end if
        call append_chunk(timetrace_luns(m),data,xy_fmt,ascii_fmt)
      else
        call output_array(trim(luname),&
           & 'diagnostic/diagnos_cross_phase', &
           & data, 'F', xy_fmt, ascii_fmt)
      end if
      if(.not.(do_append .and. exists)) then
        call attach_metadata_grid(trim(luname), &
           & 'diagnostic/diagnos_cross_phase', &
           & 'cross_phase_degr_grid', 'krho', ascii_fmt)
        call attach_metadata(trim(luname), &
           & 'diagnostic/diagnos_cross_phase', phys_unit_key, &
           & 'probability density', ascii_fmt)
        call attach_metadata(trim(luname), &
           & 'diagnostic/diagnos_cross_phase', description_key, &
           & 'PDF of the cross-phase of '// &
           & trim(tokens(cross_phases_parsed(1,i)))//' and '// &
           & trim(tokens(cross_phases_parsed(2,i))), &
           & ascii_fmt)
        call attach_metadata(trim(luname), &
           & 'diagnostic/diagnos_cross_phase', comments_key, &
           & 'This is the probability of the phase difference of the&
           & fluctuations at one position to be in the interval&
           & [angle-delta/2,angle+delta/2], as a function of the binormal&
           & wavenumber. The phase difference of the fluctuations A and B&
           & is computed as atan(Im(A/B)/Re(A/B)). Note that the &
           & probability density function is a density with respect to &
           & radiants, not degrees.', ascii_fmt)
        if(weight_mode == WEIGHTED_YSX) then
          call attach_metadata(trim(luname), &
             & 'diagnostic/diagnos_cross_phase', 'weighting', &
             & 'To make plots cleaner, the PDF is weighted to a sum over&
             & y, s and x.', ascii_fmt)
        else if(weight_mode == WEIGHTED_SX) then
          call attach_metadata(trim(luname), &
             & 'diagnostic/diagnos_cross_phase', 'weighting', &
             & 'To make plots cleaner, the PDF is weighted to a sum over&
             & s and x.', ascii_fmt)
        end if
      end if
    end if
  end subroutine output_pdf

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_degr_grid()
    use io, only : output_array, attach_metadata, xy_fmt
    use io, only : ascii_fmt, description_key, comments_key, phys_unit_key

    integer :: i
    
    call output_array('cross_phase_degr_grid',&
       & 'diagnostic/diagnos_cross_phase', &
       & (/ (-180 + (i-1)*360.0/(cross_phase_nclasses-1), &
       & i = 1, cross_phase_nclasses) /), &
       & 'F', xy_fmt, ascii_fmt)
    call attach_metadata('cross_phase_degr_grid', &
       & 'diagnostic/diagnos_cross_phase', phys_unit_key, &
       & 'degree', ascii_fmt)
    call attach_metadata('cross_phase_degr_grid', &
       & 'diagnostic/diagnos_cross_phase', description_key, &
       & 'fluctuation phase difference', &
       & ascii_fmt)
    call attach_metadata('cross_phase_degr_grid', &
       & 'diagnostic/diagnos_cross_phase', comments_key, &
       & 'Each class extends from value-delta/2 to value+delta/2, except the&
       & outer two classes which have just half this size.', ascii_fmt)
  end subroutine output_degr_grid

end module diagnos_cross_phase
