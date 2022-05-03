!------------------------------------------------------------------------------
!> Calculates growth rates, amplitudes, and frequencies for both linear
!> and nonlinear runs, using the sum of the perturbed fields.  
!> 
!> Uses mode amplitudes calculated by the normalise module.
!> Calculations are both per-mode and global.  In the per-mode case, 
!> connected kx modes are identified using the mode_label array.
!> The 1D outputs (as a function of krho) give only the kx = 0 mode.
!>
!> Performs growth rate convergence checks to trigger stop in linear runs.
!------------------------------------------------------------------------------
module diagnos_growth_freq
  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output, screen_output, read_last_data
  public :: calc_largestep, calc_smallstep

  public :: growth_rates, real_frequency

  !> Switch, to calculate the mode frequencies.
  !> output the frequencies of individual toroidal modes
  logical, save, public :: lfrequencies
  
  !> Switch to calculate and output the growth rate of the dominant
  !> mode (also called "global growth rate"), as well as the growth
  !> rates of individual modes.
  logical, save, public :: lgrowth_rates
  
  !> calculate and output the amplitudes of individual modes
  logical, save, public :: lamplitudes

  !> mode amplitudes
  real, save, allocatable, dimension(:,:) :: amplitudes
  !> mode amplitudes at previous time; used for growth rates
  real, save, allocatable, dimension(:,:) :: last_amplitudes
  logical, save, allocatable, dimension(:,:) :: amplitude_is_tiny
  !> mode amplitudes temporary array
  real, save, allocatable, dimension(:,:) :: amplitudes_G
  !> growth rates, per toroidal mode
  real, save, allocatable, dimension(:) :: growth_rates

  integer, parameter :: history_size = 6
  real, allocatable, save :: growth_rate_history(:,:)

  real, allocatable, save :: std_gr(:), mean_gr(:)

  !> the growth rate (or maximum growth rate) as reported in time.dat
  real, save :: growth_rate

  !> The real frequency of the mode.
  real, allocatable, save :: real_frequency(:)
  real, save :: single_real_frequency

  !> phase, last phase, for real frequency calculation
  real, save :: single_phase, single_last_phase
  real, allocatable, save :: phase(:),last_phase(:)
  complex, allocatable, save :: phase_tmp(:)

  !> time between phase calculations
  real, save :: delta_time_phase = 1.

  !> sum of the amplitudes to calculate the growth rate of connected
  !> modes, as for real_frequency
  real, allocatable, save :: sum_amp(:), sum_last_amp(:)

  !> logical unit numbers of various output files (new format/order/grouping)
  integer, save :: i_growth_rate, i_real_freq

  !> logical unit numbers for output files (in legacy format/order/grouping)
  integer, save :: i_time_legacy
  integer, save :: i_growth_rates, i_amplitudes, i_frequencies
  integer, save :: i_growth_rates_all_modes, i_frequencies_all_modes

  !> store the index of the current strongest growing mode. It's
  !> frequency is output as the 'dominant frequency' in case of
  !> normalisation per mode.
  integer, save :: current_dominating_mode_label
  
  !> 
  integer :: legacy_ncolumns

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    
    lamplitudes = .false.
    lgrowth_rates = .true.
    lfrequencies = .true.
    
  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lgrowth_rates,1)
    call mpibcast(lfrequencies,1)
    call mpibcast(lamplitudes,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use grid, only : nmod
    use control, only : normalize_per_toroidal_mode, ntime
    use general, only : gkw_warn

    if (nmod > 1 .and. normalize_per_toroidal_mode .and. .not. lgrowth_rates) then
      call gkw_warn('Normalisation per toroidal mode: Enable per-mode' // &
         & ' amplitudes, growth rates and frequencies diagnostics.')
      lamplitudes = .true.
      lgrowth_rates = .true.
      lfrequencies = .true.
    end if

    if(ntime == 1 .and. normalize_per_toroidal_mode) then
      call gkw_warn('Correct measurement of the dominant real frequency at the&
         & first timestep is not implemented for normalize_per_toroidal_mode=T.')
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : comments_key, description_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use grid, only : nmod
    use control, only : io_legacy, non_linear, spectral_radius
    use global, only : dotdat
    logical, intent(inout) :: requirements(:,:)

    ! Array so far not used, thus to keep the compiler quiet.
    if (.false.) write(*,*) requirements

    ! FIXME Rather than initialising with 0, one could initialise with
    ! the strongest mode, as growth rate is not yet known. For
    ! restarted runs, this is often also the strongest growing
    ! mode. This given a sensible value already at the first timestep
    ! after restarts.
    current_dominating_mode_label = 0

    if(non_linear) then
      legacy_ncolumns = 1
    else
      ! besides time, write also growth rate and real frequency into
      ! the same file
      legacy_ncolumns = 3
    end if

    if(root_processor) then
      ! TODO: check if last_largestep_time is correctly in all cases
      ! (e.g. after restart)

      ! open the logical unit to write the time evolution of the
      ! growth rate.
      call open_time_lu()

      ! open the logical units for the growth rates and freqs of the
      ! individual modes.
      if (lgrowth_rates) then
        if(io_legacy) then
          call open_real_lu('growth.dat', 'diagnostic/diagnos_growth_freq', &
             & (/ nmod /), ascii_fmt, i_growth_rates)
        else
          call open_real_lu('growth_rates', 'diagnostic/diagnos_growth_freq', &
             & (/ nmod /), ascii_fmt, i_growth_rates)
          call attach_metadata_grid(i_growth_rates, 'time', 'krho', ascii_fmt)
          call attach_metadata(i_growth_rates, phys_unit_key, &
             & 'v_{th,ref}/R_{ref}', ascii_fmt)
          call attach_metadata(i_growth_rates, description_key, &
             & 'Growth rates of the individual modes with kx=0 at the LFS.', &
             & ascii_fmt)
          call attach_metadata(i_growth_rates, comments_key, not_avail, ascii_fmt)
        end if

        if(spectral_radius) then
          call open_real_lu('growth_rates_all_modes', 'diagnostic/diagnos_growth_freq', &
             & shape(growth_rates), ascii_fmt, i_growth_rates_all_modes)
          call attach_metadata_grid(i_growth_rates_all_modes, 'time', 'mode_label', &
             & ascii_fmt)
          call attach_metadata(i_growth_rates_all_modes, phys_unit_key, &
             & 'v_{th,ref}/R_{ref}', ascii_fmt)
          call attach_metadata(i_growth_rates_all_modes, description_key, &
             & 'Growth rates of all individual modes, with any kx at the LFS.', &
             & ascii_fmt)
          call attach_metadata(i_growth_rates_all_modes, comments_key, &
             & 'Enumeration is according to mode_label.', ascii_fmt)

          call open_real_lu('frequencies_all_modes', 'diagnostic/diagnos_growth_freq', &
             & shape(real_frequency), ascii_fmt, i_frequencies_all_modes)
          call attach_metadata_grid(i_frequencies_all_modes, 'time', 'mode_label', &
             & ascii_fmt)
          call attach_metadata(i_frequencies_all_modes, phys_unit_key, &
             & 'v_{th,ref}/R_{ref}', ascii_fmt)
          call attach_metadata(i_frequencies_all_modes, description_key, &
             & 'Real frequencies of all individual modes, with any kx at the LFS.', &
             & ascii_fmt)
          call attach_metadata(i_frequencies_all_modes, comments_key, &
             & 'Enumeration is according to mode_label.', ascii_fmt)
        end if
      end if

      if (lfrequencies) then
        call open_real_lu(dotdat('frequencies',io_legacy), &
           & 'diagnostic/diagnos_growth_freq',&
           & (/ nmod /), ascii_fmt, i_frequencies)
        call attach_metadata_grid(i_frequencies, 'time', 'krho', ascii_fmt)
        call attach_metadata(i_frequencies, phys_unit_key, 'v_{th,ref}/R_{ref}', &
           & ascii_fmt)
        call attach_metadata(i_frequencies, description_key, &
           & 'Real frequencies of the individual modes', ascii_fmt)
        call attach_metadata(i_frequencies, comments_key, not_avail, ascii_fmt)
      end if

      if (lamplitudes) then
        call open_real_lu(dotdat('amplitudes',io_legacy), &
           & 'diagnostic/diagnos_growth_freq', &
           & (/ nmod /), ascii_fmt, i_amplitudes)
        call attach_metadata_grid(i_amplitudes, 'n_x_grid,time', 'krho', ascii_fmt)
        call attach_metadata(i_amplitudes, phys_unit_key, not_avail, &
           & ascii_fmt)
        call attach_metadata(i_amplitudes, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_amplitudes, comments_key, not_avail, ascii_fmt)
      end if
    end if

    real_frequency(:) = 0.0
    single_real_frequency = 0.0
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine open_time_lu()
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : comments_key, description_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use control, only : io_legacy, non_linear

    if(io_legacy) then
      ! open the file to write or read the time evolution of the growth rate.
      call open_real_lu('time.dat', 'diagnostic/diagnos_growth_freq', &
         & (/ legacy_ncolumns/), ascii_fmt, &
         & i_time_legacy)
    else
      if(.not.non_linear) then
        ! the time is output by diagnos_grid
        call open_real_lu('dominant_growth_rate', &
           & 'diagnostic/diagnos_growth_freq', (/1/), &
           & ascii_fmt, i_growth_rate)
        call attach_metadata_grid(i_growth_rate, 'time', ascii_fmt)
        call attach_metadata(i_growth_rate, phys_unit_key, 'v_{th,ref}/R_{ref}', &
           & ascii_fmt)
        call attach_metadata(i_growth_rate, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_growth_rate, comments_key, not_avail, ascii_fmt)
        
        call open_real_lu('dominant_real_freq', &
           & 'diagnostic/diagnos_growth_freq', (/1/), &
           & ascii_fmt, i_real_freq)
        call attach_metadata_grid(i_real_freq, 'time', ascii_fmt)
        call attach_metadata(i_real_freq, phys_unit_key, 'v_{th,ref}/R_{ref}', &
           & ascii_fmt)
        call attach_metadata(i_real_freq, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_real_freq, comments_key, not_avail, ascii_fmt)
      end if
    end if
  end subroutine open_time_lu

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine close_time_lu()
    use io,         only : close_lu, ascii_fmt
    use control,    only : io_legacy, non_linear

    if(io_legacy) then
      call close_lu(i_time_legacy, ascii_fmt)
    else
      if(.not.non_linear) then
        ! the time is output by diagnos_grid
        call close_lu(i_growth_rate, ascii_fmt)
        call close_lu(i_real_freq, ascii_fmt)
      end if
    end if

  end subroutine close_time_lu
  
  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general, only : gkw_abort
    use grid, only : nmod, nx, n_x_grid
    use mode, only : nmodes_G
    integer :: ierr, i

    allocate(std_gr(nmodes_G), stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: Could not allocate std_gr')
    allocate(mean_gr(nmodes_G), stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: Could not allocate mean_gr')
    mean_gr(:) = 0.0

    ! arrays for mode amplitudes growth rates
    !if (lamplitudes .or. lgrowth_rates) then
      allocate(amplitudes(nmod,nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: amplitudes')
      allocate(last_amplitudes(nmod,nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: last_amplitudes')
      allocate(amplitude_is_tiny(nmod,nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: amplitude_is_tiny')
      amplitude_is_tiny = .false.
      ! allocate the tmp array in any case because one is needed to
      ! reduce-sum the amplitudes
      allocate(amplitudes_G(nmod,n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: amplitudes_G')
    !end if

    ! array for growth rates
    if (lgrowth_rates) then
      allocate(growth_rates(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: growth_rates')
      allocate(sum_amp(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: sum_amp')
      allocate(sum_last_amp(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: sum_last_amp')
    end if

    if (nmodes_G == 0) call gkw_abort('allocate diagnostic after kgrid')

    ! mode frequency calculation arrays
    !if (lfrequencies) then
      allocate(real_frequency(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: real_frequency')
      allocate(phase(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: phase')
      allocate(last_phase(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: last_phase')
      allocate(phase_tmp(nmodes_G),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: phase_tmp')
    !end if

    ! allocate the array which holds the last few growth rates,
    ! for every toroidal mode, if requested
    allocate(growth_rate_history(nmodes_G, history_size), stat = ierr)
    if (ierr /= 0) call gkw_abort('could not allocate growth_rate_history')
    ! initialise the history array with something that does not trigger
    ! the convergence check
    do i = 1, history_size
      growth_rate_history(:, i) = i
    end do



  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine read_last_data()
    use control,      only : nt_complete, io_legacy
    use control,      only : time
    use io,           only : lu_exists, ascii_fmt
    use io,           only : read_last_chunk
    use mpiinterface, only : root_processor, mpibcast
    use general,      only : gkw_warn
    integer :: nlast
    real, dimension(3) :: last_chunk

    if(io_legacy) then
      if(.not. lu_exists('time.dat', 'grid', ascii_fmt)) then
        if(root_processor) then
          write(*,*) 'file or dataset time.dat not found. Do not read last value.'
          ! Otherwise the time value set in the restart
          ! module is kept.
        end if
      else
        
        if(root_processor) then
          ! Get the last time value of the previous run from the output data.
          call open_time_lu()
          call read_last_chunk(i_time_legacy, '(3(es13.5))', &
             & last_chunk(1:legacy_ncolumns), &
             & nlast, ascii_fmt)
          call close_time_lu()
          time = last_chunk(1)

          if(nt_complete /= nlast) then
            call gkw_warn('The logical unit "time(.dat)" '// &
               & 'is shorter or longer than expected.');
            write (*,*) "nt_complete = ", nt_complete
            write (*,*) "number of chunks in logical unit = ", nlast
          end if
        end if

        ! Broadcast the retrieved values to all processes.
        call mpibcast(time,1)
      end if
    end if

  end subroutine read_last_data

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()
    use grid, only : proc_subset

    if(proc_subset(0,0,1,1,1)) then
      ! re-calculate the mode amplitudes, to avoid a spike in the
      ! growthrate time trace
      call calc_amplitudes()
    end if

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output()

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(i_smallstep)
    use control, only : naverage
    integer, intent(in) :: i_smallstep

    if (i_smallstep == naverage - 1) then
      call calc_phase()
    else if (i_smallstep == naverage) then
      call calc_phase()
    end if

  end subroutine calc_smallstep

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use grid, only : proc_subset
    use control, only : stop_me, nan_stop, method
    use mpiinterface, only : mpiallreduce_or_inplace
    
    if(proc_subset(0,0,1,1,1)) then
      ! calculate the mode amplitudes
      !if (lamplitudes .or. lgrowth_rates) then
      call calc_amplitudes()
      !end if

      ! calulate the growth rate
      !if(lgrowth_rates) then
      call diagnostic_growth_rate()
      if(method /= 'EIV') then
        call check_convergence_and_stop()
      end if
      !end if

      ! (calculate the real frequency)
      !if (lfrequencies) then
      call diagnostic_real_freq()
      !end if
    end if

    call mpiallreduce_or_inplace(stop_me)
    call mpiallreduce_or_inplace(nan_stop)

  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use grid, only : nmod, n_x_grid
    use mode, only : ixzero, mode_label_G
    use mpiinterface, only : root_processor
    use control, only : spectral_radius
    integer :: imod, ix

    if (root_processor) then
      call write_time_output()

      if (lamplitudes) then
        do ix = 1, n_x_grid
          call append_chunk(i_amplitudes, &
             & amplitudes_G(:, ix), xy_fmt, ascii_fmt)
        end do
      end if

      if (lgrowth_rates) then
        call append_chunk(i_growth_rates, &
           & (/ (growth_rates(mode_label_G(imod,ixzero)),imod=1,nmod) /), &
           & xy_fmt, ascii_fmt)
        if(spectral_radius) then
          call append_chunk(i_growth_rates_all_modes, &
             & growth_rates, &
             & xy_fmt, ascii_fmt)
          call append_chunk(i_frequencies_all_modes, &
             & real_frequency, &
             & xy_fmt, ascii_fmt)
        end if
      end if

      ! The frequencies are only outupt for the kx=0 modes
      ! and do not exactly correspond to the frequencies above
      if (lfrequencies) then
        call append_chunk(i_frequencies, &
           & (/ (real_frequency(mode_label_G(imod,ixzero)),imod=1,nmod) /), &
           & xy_fmt, ascii_fmt)
      end if
    end if

  end subroutine output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine screen_output()
    use mpiinterface, only : root_processor
    use control,   only : non_linear
    use grid,      only : nmod, n_x_grid
    use mode,      only : mode_box, mode_label_G
    integer :: imod, ix

    ! Only root process shall write to terminal
    if(.not.root_processor) return

    if(mode_box) then
      if(.not. non_linear) then
50      format('Global growth rate : ',es13.5)
        write(*,50) growth_rate
60      format('Real frequency : ',es13.5)
        write(*,60) single_real_frequency
      end if
    else
      mod : do imod = 1, nmod
        x : do ix = 1, n_x_grid
          if(.not. non_linear) then
201         format('nMode ',I3,' xMode ',I3,' Growth rate ',es13.5)
            write(*,201)imod,ix,growth_rate
202         format('nMode ',I3,' xMode ',I3,' Real frequency ',es13.5)
            write(*,202)imod,ix,real_frequency(mode_label_G(imod,ix))
          end if
        end do x
      end do mod
    end if

  end subroutine screen_output


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Writes the given values to time.dat (legacy setting) or to three
  !> separate logical units.
  !> - For non_linear runs only the time is written.
  !> - For linear runs the frequency and the growth rate is written
  !>
  !----------------------------------------------------------------------------
  subroutine write_time_output()
    use control, only : time, non_linear, io_legacy
    use io, only : append_chunk, ascii_fmt
    use mpiinterface, only : root_processor
    if(.not. root_processor) return

    if(io_legacy) then
      if (non_linear) then
        ! For non_linear runs write time only.
        call append_chunk(i_time_legacy, (/ time /), '(3(es13.5))', ascii_fmt)
      else
        ! write time and growth rate and real frequency
        call append_chunk(i_time_legacy, &
           & (/ time, growth_rate, single_real_frequency /), &
           & '(3(es13.5))', ascii_fmt)
      end if
    else
      if(.not. non_linear) then
        ! new style: produce three separate datasets.
        call append_chunk(i_growth_rate, (/ growth_rate /), '(1(es13.5))', &
           & ascii_fmt)
        call append_chunk(i_real_freq, (/ single_real_frequency /), &
           & '(1(es13.5))', ascii_fmt)
      end if
    end if

  end subroutine write_time_output

  !---------------------------------------------------------------------------
  !> construct the array for quick look up of the fields
  !---------------------------------------------------------------------------
  function get_field_index_lookup_array() result(lindx)
    use grid, only : nx, nmod, ns
    use dist, only : iphi, iapar, ibpar, number_of_fields
    use control, only : nlphi, nlapar, nlbpar
    use index_function, only : indx
    integer, dimension(number_of_fields*ns*nx*nmod) :: lindx
    integer, allocatable, dimension(:) :: field_id
    integer :: nfields, k, ix, imod, is, idx

    ! construct the array for quick look up of the fields
    allocate(field_id(number_of_fields))
    nfields = 0
    if (nlphi) then
      nfields = nfields + 1
      field_id(nfields) = iphi
    end if
    if (nlapar) then
      nfields = nfields + 1
      field_id(nfields) = iapar
    end if
    if (nlbpar) then
      nfields = nfields + 1
      field_id(nfields) = ibpar
    end if

    idx = 0
    do k = 1, number_of_fields
      do ix = 1, nx
        do imod = 1, nmod
          do is = 1, ns
            idx=idx+1
            lindx(idx) = indx(field_id(k),imod,ix,is)
          end do
        end do
      end do
    end do
    deallocate(field_id)

  end function get_field_index_lookup_array

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine calc_amplitudes
    use mpicomms,       only : COMM_S_NE, COMM_S_NE_X_NE
    use mpiinterface,   only : mpireduce_sum_inplace
    use mpiinterface,   only : mpiallreduce_sum_inplace
    use grid,           only : ns,nmod,nx
    use dist,           only : fdis => fdisi, number_of_fields
    use control,        only : spectral_radius
    use mode,           only : iyzero

    logical, save :: first_time = .true.

    integer :: idx,i,imod,ix, ifield

    integer, save, allocatable, dimension(:) :: lindx
    real, save, allocatable, dimension(:) :: cphi

    if (first_time) then
      allocate(lindx(number_of_fields*ns*nx*nmod))
      lindx = get_field_index_lookup_array()

      if (.not. spectral_radius) allocate(cphi(nmod))

      ! initialise the amplitudes
      amplitudes(:,:) = 1.

      ! only run this piece of code once
      first_time = .false.

    end if

    last_amplitudes = amplitudes

    ! calculate the local amplitudes
    amplitudes = 0.
    idx = 0

    if (spectral_radius) then
      do ifield = 1, number_of_fields
        do ix = 1, nx
          do imod = 1, nmod
            ! ignore contributions of the zero mode
              do i = 1, ns
                idx=idx+1
                if (.not. imod == iyzero) then
                  amplitudes(imod,ix) = amplitudes(imod,ix) + abs(fdis(lindx(idx)))**2
                end if
              end do
            end do
        end do
      end do

      call mpireduce_sum_inplace(amplitudes, &
         & shape(amplitudes), COMM_S_NE)
      amplitudes = sqrt(amplitudes)
      amplitudes_G = amplitudes

    else
      cphi(:)=0.
      ! normalise only with the electrostatic potential
      do ix = 1, nx
        do imod = 1, nmod
          do i = 1, ns
            idx=idx+1
            cphi(imod) = cphi(imod) + abs(fdis(lindx(idx)))**2
          end do
        end do
      end do
      call mpiallreduce_sum_inplace(cphi, size(cphi,1),COMM_S_NE_X_NE)
      
      do imod = 1, nmod
        ! all x values are set the same for growth_rates
        amplitudes(imod,:) =  sqrt(cphi(imod))
        amplitudes_G(imod,:) =  sqrt(cphi(imod))
      end do
    end if

  end subroutine calc_amplitudes

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_phase

    use mpiinterface,   only : mpiallreduce_sum_inplace
    use mpicomms,       only : COMM_S_NE_X_NE
    use constants,      only : pi
    use control,        only : time, normalize_per_toroidal_mode
    use grid,           only : parallel_s, parallel_x, ns, nmod, nx, gx
    use dist,           only : fdisi,number_of_fields,iphi,iapar,ibpar
    use mode,           only : mode_label_G,nmodes_G
    use index_function, only : indx
    use general,        only : gkw_abort

    ! keep the last time this subroutine was called
    real, save :: last_time_diagn_phase = 0.
    integer, allocatable, save :: lindx(:,:,:)
    integer :: idx,i,j,k,is
    logical, save :: first_call = .true.
    integer :: ierr
    !> a helper array, used to iterate over fields
    integer, dimension(3) :: field_id = (/ iphi, iapar, ibpar /)

    if (first_call) then

      last_phase = 0.
      single_last_phase = 0.
      phase = 0.
      single_phase = 0.
      single_real_frequency = 0.
      phase_tmp = 0.
      allocate(lindx(nmod,nx,number_of_fields*ns),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: lindx')

      ! fill the array
      idx = 0
      lindx(:,:,:) = 0
      do k = 1, number_of_fields ; do is = 1, ns
        idx = idx + 1
        do j = 1, nx ; do i = 1, nmod
          lindx(i,j,idx) = indx(field_id(k),i,j,is)
        end do; end do
      end do; end do

      first_call = .false.

    end if

    ! save time interval and update last_time_diagn_phase
    delta_time_phase = time - last_time_diagn_phase
    last_time_diagn_phase = time

    ! before calculating the new phase, store the old one
    last_phase  = phase
    single_last_phase = single_phase
    phase_tmp = (0.,0.)

    ! All regular fields are used, calculation.
    do k = 1, number_of_fields*ns
      do j = 1, nx
        do i = 1, nmod
          phase_tmp(mode_label_G(i,gx(j))) = phase_tmp(mode_label_G(i,gx(j))) + &
             & fdisi(lindx(i,j,k))
        end do
      end do
    end do

    ! sum-reduce the over the s-and x direction
    if (parallel_s .or. parallel_x) then
      call mpiallreduce_sum_inplace(phase_tmp,nmodes_G,COMM_S_NE_X_NE)
    end if

    ! implicit loop over modes
    phase = atan2(aimag(phase_tmp),real(phase_tmp))

    if(normalize_per_toroidal_mode) then
      ! either take the phase of the dominant mode
      if(current_dominating_mode_label > 0) then
        single_phase = atan2(aimag(phase_tmp(current_dominating_mode_label)),&
           & real(phase_tmp(current_dominating_mode_label)))
      else
        single_phase = 0
      end if

      ! or try to reproduce the ordinary computation of the dominant
      ! phase, sort of, but I did not succeed in this and it is
      ! converges slower anyway.
    else
      single_phase = atan2(aimag(sum(phase_tmp)),real(sum(phase_tmp)))
    end if

    ! Treat the cases of 2\pi jumps.
    ! implicit loop over modes
    where (abs(phase-last_phase) > max(pi/4., &
       & abs(3.*real_frequency*delta_time_phase)))

      last_phase = last_phase - sign(2*pi,last_phase)
    end where
    if (abs(single_phase-single_last_phase) > max(pi/4., &
       & abs(3.*single_real_frequency*delta_time_phase))) then

      single_last_phase = single_last_phase - sign(2*pi,single_last_phase)
    end if


  end subroutine calc_phase

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diagnostic_real_freq

    use global,       only : r_tiny
    use control,      only : time, last_largestep_time

    ! last_largestep_time is used for the calculation of the real frequency.

    if (abs(time-last_largestep_time) < r_tiny) then
      real_frequency = 0.
      single_real_frequency = 0.
    else
      ! implicit loop over all modes
      real_frequency = (phase - last_phase) / delta_time_phase
      single_real_frequency = (single_phase - single_last_phase) &
         & / delta_time_phase
    end if

  end subroutine diagnostic_real_freq

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine finalize
    use io, only : close_lu, ascii_fmt
    use control, only : spectral_radius

    if (allocated(std_gr)) deallocate(std_gr)
    if (allocated(mean_gr)) deallocate(mean_gr)

    if (allocated(amplitudes)) deallocate(amplitudes)
    if (allocated(last_amplitudes)) deallocate(last_amplitudes)
    if (allocated(amplitude_is_tiny)) deallocate(amplitude_is_tiny)
    if (allocated(amplitudes_G)) deallocate(amplitudes_G)

    if(lgrowth_rates) then
      if (allocated(growth_rates)) deallocate(growth_rates)
      if (allocated(sum_amp)) deallocate(sum_amp)
      if (allocated(sum_last_amp)) deallocate(sum_last_amp)
    end if

    if (allocated(real_frequency)) deallocate(real_frequency)
    if (allocated(phase)) deallocate(phase)
    if (allocated(last_phase)) deallocate(last_phase)
    if (allocated(phase_tmp)) deallocate(phase_tmp)
    
    ! deallocate the private arrays and
    ! do nothing else.
    if (allocated(growth_rate_history)) deallocate(growth_rate_history)

    if(lgrowth_rates) then
      call close_lu(i_growth_rates, ascii_fmt)
      if(spectral_radius) then
        call close_lu(i_growth_rates_all_modes, ascii_fmt)
        call close_lu(i_frequencies_all_modes, ascii_fmt)
      end if
    end if
    if(lamplitudes) call close_lu(i_amplitudes, ascii_fmt)
    if(lfrequencies) call close_lu(i_frequencies, ascii_fmt)
    call close_time_lu
 
  end subroutine finalize

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine diagnostic_growth_rate()
    use global,       only : r_tiny
    use control,      only : normalized, normalize_per_toroidal_mode
    use normalise,    only : amp, amp_per_mode
    use control,      only : time, last_largestep_time, spectral_radius
    use grid,         only : nx, nmod, gx
    use mode,         only : mode_label_G, nmodes_G, iyzero
    use mpicomms, only : COMM_X_NE
    use mpiinterface, only : mpiallreduce_sum_inplace
    use general, only : gkw_warn
    use global, only : int2char
    use normalise, only : nonlin_norm_fac
    use control, only : non_linear

    integer :: imod, ix, i, j, mlabel
    integer, save :: igrowth=0
    real, save :: amp_last=1.0
    real :: one_over_nonlin_norm_fac=1.0
    
    if(non_linear) then
      one_over_nonlin_norm_fac = 1./nonlin_norm_fac
    endif

    ! calculate the individual growth rates if requested
    if (lgrowth_rates) then
      sum_amp = 0.
      sum_last_amp = 0.

      if (abs(time-last_largestep_time) < r_tiny) then
        ! Assume zero growth for small time changes.
        growth_rates(:) = 0.
      else 
        ! note that in nonspectral case every x has same mode label and
        ! same amplitude
        do ix=1,nx
          do imod=1,nmod
            mlabel = mode_label_G(imod,gx(ix))
            sum_amp(mlabel) = sum_amp(mlabel) + amplitudes(imod,ix)**2
            sum_last_amp(mlabel) = sum_last_amp(mlabel) &
               & + last_amplitudes(imod,ix)**2
          end do
        end do

        ! reduce to *all*, to be able to check convergence on any process
        call mpiallreduce_sum_inplace(sum_amp, nmodes_G, COMM_X_NE)
        call mpiallreduce_sum_inplace(sum_last_amp, nmodes_G, COMM_X_NE)
        
        do ix=1,nx
          do imod=1,nmod
            mlabel = mode_label_G(imod,gx(ix))
            ! if the mode is not vanishingly small
            if(abs(sum_amp(mlabel)) >= r_tiny .and. abs(sum_last_amp(mlabel)) >= r_tiny) then
              if(normalized) then
                growth_rates(mlabel) = log(amp*sqrt(sum_amp(mlabel)/sum_last_amp(mlabel))) / &
                   & (time-last_largestep_time)
              else if(normalize_per_toroidal_mode) then
                growth_rates(mlabel) = log(amp_per_mode(imod)*one_over_nonlin_norm_fac*&
                   & sqrt(sum_amp(mlabel)/sum_last_amp(mlabel))) / &
                   & (time-last_largestep_time)
              else
                ! not normalised - dont need to use the factors
                growth_rates(mlabel) = log(sqrt(sum_amp(mlabel)/sum_last_amp(mlabel))) / &
                   & (time-last_largestep_time)
              end if
            else
              if(.not. amplitude_is_tiny(imod,ix)) then
                if(imod /= iyzero) then
                  ! log a warning when normalisation quenches modes:
                  ! spit out this information only a single time, and
                  ! never for the zero-mode (because it is not
                  ! linearly unstable anyway)
                  call gkw_warn('amplitude of mode imod='//trim(int2char(imod,3))//&
                     & ' ix='//trim(int2char(ix,3))//' mode_label='//trim(int2char(mlabel,4))//&
                     & ' has become tiny, cannot measure its growthrate anymore.')
                end if
                amplitude_is_tiny(imod,ix) = .true.
              end if
              growth_rates(mlabel) = 0
            end if
          end do
        end do

      end if

      current_dominating_mode_label = maxloc(growth_rates,1)
      
    end if

    ! calculate the dominant growth rate which is output to time.dat by
    ! default.  Nota Bene need to use amp here.  The dominant growth_rate is
    ! calculated regardless of whether individual growth rates are
    ! requested with lgrowth_rates.
    
    if (abs(time-last_largestep_time) < r_tiny) then
      growth_rate = 0.
    else
      if(normalize_per_toroidal_mode) then
        growth_rate = maxval(growth_rates)
      else
        growth_rate = (log(amp)-log(amp_last)) / (time - last_largestep_time)
        if(.not.normalized) then
          amp_last = amp
        end if
      end if
    end if

    ! In what follows, calculate quantities used for the convergence check:

    ! fill in the array storing the growth rate time history
    igrowth = mod(igrowth+1,history_size)
    if (lgrowth_rates .and. spectral_radius) then
      growth_rate_history(:, igrowth+1) = growth_rates
    else
      ! FIXME in the nonspectral case the value of growth_rate and
      ! growth_rates differ (cf. testcase input_non_spectral.dat).
      ! As they are not consistent in that case at the moment, fall
      ! back to the global growth_rate for the convergendce check,
      ! although there are actually nmodes_G > 1 individual modes in
      ! and lgrowth_rates=.true. .
      growth_rate_history(:, igrowth+1) = growth_rate
    end if

    ! standard deviation of the growth rate history
    std_gr = 0.
    do i = 1, nmodes_G
      mean_gr(i) = sum(growth_rate_history(i,:)) / history_size
    end do

    do i = 1, nmodes_G
      do j = 1, history_size
        std_gr(i) = std_gr(i) + (growth_rate_history(i,j)-mean_gr(i))**2
      end do
      std_gr(i) = sqrt(std_gr(i) / history_size)
    end do

  end subroutine diagnostic_growth_rate

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine check_convergence_and_stop
    use control, only : gamatol, stop_me, nan_stop, min_gr
    use control, only : max_gr, neoclassics
    use control, only : time, normalize_per_toroidal_mode
    use mpiinterface, only : root_processor
    use mode, only : nmodes_G
    integer :: i
    
    if (gamatol > 0. .and. .not. stop_me) then
      if(normalize_per_toroidal_mode) then
        ! check for convergence of the growth rates of *only the
        ! unstable modes*

        ! if there is at least one stable mode
        stop_me = count(mean_gr > 0) > 0

        do i = 1, nmodes_G
          if(mean_gr(i) > 0) then
            stop_me = stop_me .and. (std_gr(i) < gamatol)
          end if
        end do
      else
        ! check for convergence of the growth rates of *all modes*
        stop_me = (maxval(std_gr) < gamatol) .or. stop_me
      end if
      if (stop_me .and. root_processor) then
        if(lgrowth_rates) then
          if (normalize_per_toroidal_mode) then
            write(*,*) '**** Growth rate convergence reached for all growing modes: STOP ****'
          else
            write(*,*) '**** Growth rate convergence reached for all modes: STOP ****'
          end if
        else
          write(*,*) '**** Growth rate convergence reached: STOP ****'
        end if
      end if
    end if
      
    if (maxval(mean_gr) > max_gr .and. time > 1.0 .and. .not. stop_me) then
      if(root_processor) then
        write(*,*) '**** max_gr: Run unstable: STOP ****'
        write(*,*)
      end if
      stop_me = .true.
      nan_stop = .true.
    end if

    if (maxval(mean_gr) < min_gr .and. maxval(std_gr) < 10.*gamatol &
       & .and. .not. neoclassics &
       & .and. .not. stop_me) then
      if(root_processor) then
        if(lgrowth_rates) then
          write(*,*) '**** min_gr: All modes are stable: STOP ****'
        else
          write(*,*) '**** min_gr: Mode is stable: STOP ****'
        end if
        write(*,*)
      end if
      stop_me = .true.
    end if

  end subroutine check_convergence_and_stop

end module diagnos_growth_freq
