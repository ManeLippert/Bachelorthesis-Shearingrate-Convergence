!------------------------------------------------------------------------------
!> Outputs selected points and xy slices of the distribution function.
!------------------------------------------------------------------------------
module diagnos_f

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  logical, save, public :: loutput_fdis
  
  !> Time interval of temporal average of the distribution
  integer, save, public :: tavg_start, tavg_end, tavg_len

  !> Private FFT arrays.
  complex, save, allocatable :: a(:,:)
  real, save, allocatable :: ar(:,:)

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
  
    use control, only : ntime
    
    loutput_fdis = .true.
    
    ! average over whole simulation
    tavg_start = 1
    tavg_end = ntime
    
  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(loutput_fdis,1)
    call mpibcast(tavg_start,1)
    call mpibcast(tavg_end,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : ntime
    use general, only : gkw_warn
     
    if(.not.loutput_fdis) return
    
    if(tavg_start < 1) then
      call gkw_warn('For temporal average of the distribution tavg_start < 1 makes no sense. &
                    & tavg_start = 1 is set! ')
      tavg_start = 1
    end if
    
    if(tavg_start > ntime) then
      call gkw_warn('For temporal average of the distribution tavg_start > makes no sense. &
                    & tavg_start = ntime is set! ')
      tavg_end = ntime
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface,    only : register_tag_range
    use grid,            only : number_of_species
    use control,         only : laverage_dist_over_time
    use global,          only : DISTRIBUTION
    use diagnos_generic, only : LOCAL_DATA
    
    logical, intent(inout) :: requirements(:,:)
    
    if(.not.loutput_fdis .and. .not. laverage_dist_over_time) return
    
    if(laverage_dist_over_time) then
    
      !> set requirements
      requirements(DISTRIBUTION,LOCAL_DATA) = .true.
      
      !> length of averaging time interval
      tavg_len = tavg_end-tavg_start+1
    
    end if
    

    call register_tag_range(number_of_species*2, &
       & tag_range_start, tag_range_end_inkl)
       
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general,          only : gkw_abort
    use diagnos_generic, only : mphiw3t,mphit,mrad_l
    integer :: ierr
    
    if(.not.loutput_fdis) return

    !Private FFT arrays for diagnostics only.
    allocate(a(mphiw3t,mrad_l), stat = ierr) 
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate a in diagnostic')
    endif

    allocate(ar(mphit,mrad_l), stat = ierr) 
    if (ierr.ne.0) then 
      call gkw_abort('Could not allocate ar in diagnostic')
    endif

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    if(.not.loutput_fdis) return

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()
    if(.not.loutput_fdis) return

  end subroutine initial_output

  !--------------------------------------------------------------------

  !>
  !--------------------------------------------------------------------
  subroutine final_output(number)
    integer, optional, intent(in) :: number
    
    if(.not.loutput_fdis) return

    ! For the case of multiple eigenmodes
    if (present(number)) then
      ! to number parallel.dat and velocity outputs
      call output_distr_at_point(number)
    else
      call output_distr_at_point()
    endif

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    use control, only : laverage_dist_over_time
    
    integer, intent(in) :: file_count

    ! keep compiler quiet
    if (file_count > 0) continue

    if(.not.loutput_fdis .or. .not. laverage_dist_over_time) return

    !call fdisi_xy_output(n_s_grid/2, n_mu_grid/2, n_vpar_grid/2, 1, file_count)
    
    if(laverage_dist_over_time) then
      call average_dist_over_time()
    end if

  end subroutine output

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Output fdisi x-y slice at a given mu, vpar, species point
  !> both spectral and non-spectral
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine fdisi_xy_output(s_point, mu_point, vpar_point, &
     & species, file_count)

    use general,          only : gkw_abort
    use global,           only : int2char_zeros
    use grid,             only : nmod, nx, ns, n_x_grid
    use dist,             only : fdisi, ifdis, phi
    use index_function,   only : indx
    use non_linear_terms, only : jind, mphi, nl_initialised
    use mpicomms,         only : COMM_X_NE, COMM_DUMMY
    use mpiinterface,     only : gather_array
    use diagnos_generic,  only : four2real_2D, mrad_G, mrad_l
    use fft, only : FFT_INVERSE

    integer, intent(in) :: mu_point,vpar_point,species, s_point
    integer, intent(in) :: file_count

    real :: rdum(mphi,mrad_G)
    complex :: cdum(mphi,n_x_grid)
    integer :: i,j,imod,ix,ipar
    character (len=16) :: luname

    ! Check if the initialization has been done 
    if (.not.nl_initialised) then 
      call gkw_abort ('xy_output: First call nonlinear_init')
    end if

    ! phi is used as a dummary array
    do i=1,ns; do imod = 1, nmod; do ix = 1, nx
      phi(imod,ix,i)=fdisi(indx(ifdis,imod,ix,i,mu_point,vpar_point,species)) 
    end do; end do; end do

    ! Select the point on the field line for the plot  
    do ipar = s_point, s_point 

      ! gather global array in x if needed
      call gather_array(cdum(1:nmod,1:n_x_grid),nmod,n_x_grid,  &
         & phi(1:nmod,1:nx,ipar),nmod,nx,        &
         & COMM_DUMMY,COMM_X_NE,ALLGATHER=.true.)

      luname = 'spc_fdis'//int2char_zeros(file_count,6)
      open(13,file = luname) 
      do imod = 1, nmod 
        write(13,11)(abs(cdum(imod,ix)),ix = 1, n_x_grid) 
      end do
      close(13)

      a = (0.,0.)
      ar = 0.

      do imod = 1, nmod; do ix = 1, nx
        a(imod,jind(ix)) = phi(imod,ix,ipar)
      end do; end do   

      ! Do the inverse FFT (also for non spectral)
      call four2real_2D(ar,a,FFT_INVERSE)

      ! gather global array in x if needed
      call gather_array(rdum(1:mphi,1:mrad_G),mphi,mrad_G,  &
         & ar(1:mphi,1:mrad_l),mphi,mrad_l,  &
         & COMM_DUMMY,COMM_X_NE,ALLGATHER=.true.)

      ! temporary output of fdis 
      luname = 'fdis'//int2char_zeros(file_count,6)
      open(13,file=luname) 
      do i = 1, mphi 
        write(13,11)(rdum(i,j), j = 1, mrad_G)
11      format(1024(es13.5,1x))
      end do
      close(13) 

    end do

    return 
  end subroutine fdisi_xy_output
  
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Subroutine that averages fdisi over a given interval of large time
  !> steps [tavg_start, tavg_end] of length tavg_len.
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine average_dist_over_time()
  
    use control, only : itime
    use dist,    only : fdisi, fdisi_tavg, nsolc
    use restart, only : write_restart_file
    
    integer :: i
    
    ! check if current time step lies in averaging interval
    if(itime >= tavg_start .and. itime <= tavg_end) then
    
      ! do averaging (fdisi_tavg has been initialized to zero beforehand)
      do i = 1, nsolc
        fdisi_tavg(i) = fdisi_tavg(i) + fdisi(i)/tavg_len
      end do
    
    end if
     
  end subroutine average_dist_over_time
  

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Output velocity space diagnostic for all species at s point xy_slice_ipar.
  !> Assumes grids already written by velocity_space_output.
  !> At present only useful for mode_box false runs - plots mode(1,1).
  !--------------------------------------------------------------------
  subroutine output_distr_at_point(number)

    use global,         only : int2char_zeros
    use dist,           only : ifdis, fdisi
    use index_function, only : indx
    use grid,           only : nmu,nvpar, ls, lsp, lrx, number_of_species
    use grid,           only : nmod, proc_subset
    use velocitygrid,   only : intmu,intvp
    use mode,           only : ixzero
    use diagnos_generic, only : velocity_slice_output, xy_slice_ipar
    use io, only : binary_fmt

    integer, optional, intent(in) :: number

    integer :: i, j, isl, ispg, ispl, ixg
    real, allocatable, dimension(:,:) :: local_vpar_mu
    character (len=31) ::  luname

    if (nmod > 1) return  !only output fdisi for linear runs, single mode.

    ! set ix here to ixzero 
    ixg = ixzero 

    ! allocate arrays to contain the full slice and local slice
    allocate(local_vpar_mu(nvpar,nmu))

    !USE xy_slice_ipar as for xy-slices
    !local s
    isl=ls(xy_slice_ipar)

    !Loop over all species
    species: do ispg=1,number_of_species

      if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then

        !local species
        ispl=lsp(ispg)

        ! local processor imaginary part of fdisi
        do j=1,nmu
          do i=1,nvpar
            local_vpar_mu(i,j)=aimag(fdisi(indx(ifdis,1,lrx(ixg),isl,j,i,ispl)) &
               & *intmu(j)*intvp(isl,j,i,ispl))
          end do
        end do
      end if

      ! number multiple eigenmodes if necessary
      luname="distr_at_point.imag.sp"//trim(int2char_zeros(ispg,2))
      if (present(number)) then
        luname = trim(luname)//"_eim"//trim(int2char_zeros(number,3))
      end if

      call velocity_slice_output('diagnostic/diagnos_f', &
         & local_vpar_mu,luname,ixg,xy_slice_ipar,ispg, binary_fmt, &
         & tag_range_start - 1 + (ispg-1)*number_of_species + 1)

      if(proc_subset(ixg,xy_slice_ipar,0,0,ispg)) then
        ! local processor real part of fdisi
        do j=1,nmu
          do i=1,nvpar
            local_vpar_mu(i,j)=real(fdisi(indx(ifdis,1,lrx(ixg),isl,j,i,ispl)) &
               & *intmu(j)*intvp(isl,j,i,ispl))
          end do
        end do
      end if
      ! number multiple eigenmodes if necessary
      luname="distr_at_point.real.sp"//trim(int2char_zeros(ispg,2))
      if (present(number)) then
        luname = trim(luname)//"_eim"//trim(int2char_zeros(number,3))
      end if

      call velocity_slice_output('diagnostic/diagnos_f', &
         & local_vpar_mu,luname,ixg,xy_slice_ipar,ispg, binary_fmt, &
         & tag_range_start - 1 + (ispg-1)*number_of_species + 2)

    end do species

    ! deallocate the temporary arrays
    if (allocated(local_vpar_mu)) deallocate(local_vpar_mu)

  end subroutine output_distr_at_point

end module diagnos_f
