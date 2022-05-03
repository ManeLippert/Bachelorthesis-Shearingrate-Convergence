!------------------------------------------------------------------------------
!> This diagnostic calculates the shear rate of the zonal flows.
!!    \f$ \omega = \frac{\partial}{\partial \psi} v_{E,\zeta}
!!          = \frac{\partial^2}{\partial \psi^2} \ga{\phi} \f$
!!
!! It is used for example in
!!    * Dannert and Jenko, PoP 12, 072309 (2005)
!!    * Merz and Jenko, Nucl. Fusion 50, 054005 (2010)
!!
!! Other interesting works in this context include
!!    * Hahm et al., PoP 6, No. 3, 922 (1999)
!!
!<------------------------------------------------------------------------------
module diagnos_zfshear

  implicit none

  private

  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output, final_output
  public :: output

  !> The general on/off switch for this diagnostic
  logical, save, public :: lcalc_zfshear

  !> The logical unit numbers, used to output data
  integer, save, allocatable :: i_zfshear_ga(:), i_zfflow_ga(:), i_zfphi_ga(:)
  integer :: i_zfshear, i_zfflow, i_zfphi

  !> Help array for first and second radial derivative of the electrostatic
  !> potential; the corresponding functions (dfielddx, d2fielddx2) only
  !> accept complex arrays, although the potential is a real quantity.
  !> Therefore, in the calculation of the individual quantities below, dum
  !> is casted to a real value always.
  complex, allocatable, dimension(:) :: dum
  
  !> Help arrays for the fluxsurface average as well as velocity space integral
  real, allocatable, dimension(:) :: local_x_line_ga, local_x_line
  
  !> Global arrays in radial direction containing the gyro-averaged quantities
  real, allocatable, dimension(:) :: zfshear_ga_G, zfflow_ga_G, zfphi_ga_G
  
  !> Global arrays in radial direction containing the non-gyro-averaged quantities
  real, allocatable, dimension(:) :: zfshear_G, zfflow_G, zfphi_G
  
  logical, parameter :: test_x_derivatives = .false.
  
  !> mpi datatypes 
  integer, save :: mpi_dtype_radial_line

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    ! Usually the diagnostic should be switched
    ! off by default.

    lcalc_zfshear = .false.
  end subroutine set_default_nml_values

  
  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_zfshear,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use mode, only : iyzero
    use general, only : gkw_warn

    if (.not.lcalc_zfshear) return

    if (iyzero <= 0) then
      call gkw_warn('Cannot use the zfshear diagnostic: There is no k_y=0 mode.')
      lcalc_zfshear=.false.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor, MPIREAL_X
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid, LOCAL_DATA, X_GHOSTCELLS
    use global, only : PHI_GA_FIELD, PHI_FIELD, id_x
    use control, only : spectral_radius
    use grid, only : number_of_species, proc_subset
    use global,       only : int2char_zeros
    use mpidatatypes, only : create_subarray_datatype
    
    logical, intent(inout) :: requirements(:,:)
    integer :: is
    character (len=18) :: filename

    if (.not.lcalc_zfshear) return

    requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    requirements(PHI_FIELD,LOCAL_DATA) = .true.
  
    ! requirements for radial derivatives in combination with 
    ! non-spectral runs
    if(.not.spectral_radius) then
      requirements(PHI_GA_FIELD,X_GHOSTCELLS) = .true.
      requirements(PHI_FIELD,X_GHOSTCELLS) = .true.
    end if
    
    
    ! open logical units of files containing gyro-averaged quantities
    do is = 1, number_of_species
      if (root_processor) then
      
      filename="zfshear_ga"//trim(int2char_zeros(is,2))
      call open_real_lu(trim(filename), 'diagnostic/diagnos_zfshear', &
         & shape(zfshear_ga_G), ascii_fmt, i_zfshear_ga(is))
      call attach_metadata_grid(i_zfshear_ga(is), 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfshear_ga(is), phys_unit_key, 'v_{th,ref}/R_{ref}', ascii_fmt)
      call attach_metadata(i_zfshear_ga(is), description_key,'This file has &
      & n_x_grid columns and contains the zonal flow shearing rate, i. e., the &
      & second radial derivative of the zonal part of electrostatic potential. &
      & The gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfshear_ga(is), comments_key, not_avail, ascii_fmt)

      filename="zfflow_ga"//trim(int2char_zeros(is,2))
      call open_real_lu(trim(filename), 'diagnostic/diagnos_zfshear', &
         & shape(zfshear_ga_G), ascii_fmt, i_zfflow_ga(is))
      call attach_metadata_grid(i_zfflow_ga(is), 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfflow_ga(is), phys_unit_key, '\rho_\ast v_{th,ref}', ascii_fmt)
      call attach_metadata(i_zfflow_ga(is), description_key,'This file has &
      & n_x_grid columns and contains the zonal flow velocity, i. e., the &
      & first radial derivative of the zonal part of electrostatic potential. &
      & The gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfflow_ga(is), comments_key, not_avail, ascii_fmt)

      filename="zfphi_ga"//trim(int2char_zeros(is,2))
      call open_real_lu(trim(filename), 'diagnostic/diagnos_zfshear', &
         & shape(zfshear_ga_G), ascii_fmt, i_zfphi_ga(is))
      call attach_metadata_grid(i_zfphi_ga(is), 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfphi_ga(is), phys_unit_key, '\rho_\ast^2 T_{ref}/e', ascii_fmt)
      call attach_metadata(i_zfphi_ga(is), description_key, 'This file has &
      & n_x_grid columns and contains the zonal part of the electorstatic potential. &
      & The gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfphi_ga(is), comments_key, not_avail, ascii_fmt)
      end if
    end do
      
    ! open logical units of files containing non-gyro-averaged quantities
    if(root_processor) then
      call open_real_lu('zfshear', 'diagnostic/diagnos_zfshear', &
         & shape(zfshear_G), ascii_fmt, i_zfshear)
      call attach_metadata_grid(i_zfshear, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfshear, phys_unit_key, 'v_{th,ref}/R_{ref}', ascii_fmt)
      call attach_metadata(i_zfshear, description_key, 'This file has &
      & n_x_grid columns and contains the zonal flow shearing rate, i. e., the &
      & second radial derivative of the zonal part of electrostatic potential. &
      & The non-gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfshear, comments_key, not_avail, ascii_fmt)

      call open_real_lu('zfflow', 'diagnostic/diagnos_zfshear', &
         & shape(zfflow_G), ascii_fmt, i_zfflow)
      call attach_metadata_grid(i_zfflow, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfflow, '\rho_\ast v_{th,ref}', not_avail, ascii_fmt)
      call attach_metadata(i_zfflow, description_key, 'This file has &
      & n_x_grid columns and contains the zonal flow velocity, i. e., the &
      & first radial derivative of the zonal part of electrostatic potential. &
      & The non-gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfflow, comments_key, not_avail, ascii_fmt)

      call open_real_lu('zfphi', 'diagnostic/diagnos_zfshear', &
         & shape(zfphi_G), ascii_fmt, i_zfphi)
      call attach_metadata_grid(i_zfphi, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zfphi, phys_unit_key, '\rho_\ast^2 T_{ref}/e', ascii_fmt)
      call attach_metadata(i_zfphi, description_key, 'This file has &
      & n_x_grid columns and contains the zonal part of the electorstatic potential. &
      & The non-gyro-averaged electrostatic potential is used.', ascii_fmt)
      call attach_metadata(i_zfphi, comments_key, not_avail, ascii_fmt)
      
    end if
    
    if (proc_subset(0,1,1,1,0)) then
      ! create subarray MPI datatype, used for gathering the radial
      ! profiles to root
      call create_subarray_datatype(MPIREAL_X, mpi_dtype_radial_line, &
         & id_x)
    end if
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general, only : gkw_abort
    use grid, only : nx, n_x_grid, number_of_species, proc_subset
    integer :: ierr

    if (.not.lcalc_zfshear) return
    
    
    allocate(i_zfshear_ga(number_of_species),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate i_zfshear_ga')
    allocate(i_zfflow_ga(number_of_species),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate i_zfshear_ga')
    allocate(i_zfphi_ga(number_of_species),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate i_zfshear_ga')

    allocate(dum(nx),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate dum')
    allocate(local_x_line(nx),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate local_x_line')
    allocate(local_x_line_ga(nx),stat=ierr)
    if (ierr /= 0) &
       & call gkw_abort('diagnos_zfshear :: could not allocate local_x_line_ga')
    
    
    if(proc_subset(0,1,1,1,0)) then
      allocate(zfshear_ga_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
      allocate(zfflow_ga_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
      allocate(zfphi_ga_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
    
    end if
    
    if(proc_subset(0,1,1,1,1)) then
      allocate(zfshear_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
      allocate(zfflow_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
      allocate(zfphi_G(n_x_grid),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnos_zfshear :: could not allocate')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, ascii_fmt
    use mpiinterface, only : root_processor
    use grid, only : number_of_species
    
    integer :: is

    if (.not.lcalc_zfshear) return
    
    ! deallocate all arrays of this diagnostic
    if(allocated(dum)) deallocate(dum)
    if(allocated(local_x_line)) deallocate(local_x_line)
    if(allocated(local_x_line_ga)) deallocate(local_x_line_ga)
    
    if(allocated(zfshear_ga_G)) deallocate(zfshear_ga_G)
    if(allocated(zfflow_ga_G)) deallocate(zfflow_ga_G)
    if(allocated(zfphi_ga_G)) deallocate(zfphi_ga_G)
    
    if(allocated(zfshear_G)) deallocate(zfshear_G)
    if(allocated(zfflow_G)) deallocate(zfflow_G)
    if(allocated(zfphi_G)) deallocate(zfphi_G)
    


    ! be nice and close all logical units
    do is = 1, number_of_species
      if(root_processor) then
        call close_lu(i_zfshear_ga(is), ascii_fmt)
        call close_lu(i_zfflow_ga(is), ascii_fmt)
        call close_lu(i_zfphi_ga(is), ascii_fmt)
      end if
    end do
   
    if(root_processor) then
      call close_lu(i_zfshear, ascii_fmt)
      call close_lu(i_zfflow, ascii_fmt)
      call close_lu(i_zfphi, ascii_fmt)  
    end if
  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()
  
    if (.not.lcalc_zfshear) return
    
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

    if (.not.lcalc_zfshear) return

    if (present(number)) then
      ! output for every eigenmode
    else
      ! one single output
    end if
    
  end subroutine final_output

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
    use grid, only : nx, n_x_grid, ns, nmu, nsp, number_of_species
    use grid, only : nvpar, vpmax, proc_subset, mumax, lsp
    use fields, only : get_averaged_phi
    use geom, only : ints, bn
    use velocitygrid, only : intmu, intvp
    use mode, only : iyzero
    use mpicomms, only : COMM_X_NE
    use mpicomms, only : COMM_S_NE, COMM_SP_EQ_X_EQ
    use mpiinterface, only : mpireduce_sum_inplace, gather_array
    use diagnos_generic, only : dfielddx, d2fielddx2
    use global, only : PHI_GA_FIELD, PHI_FIELD
    use dist, only : fdisi, phi, fmaxwl
    use constants, only : pi
    use io, only : ascii_fmt, xy_fmt, append_chunk
    use mpiinterface, only : root_processor

    integer :: ipar, j, is, ix, k
    real :: d3v
    

    if (.not.lcalc_zfshear) return

    ! Since gyro-averaged quantities are species dependent, calculate
    ! zfshear_ga, zfflow_ga and zfphi_ga for each species and make output.
    loop_global_species: do is = 1, number_of_species

      ! gyro-averaged potential
      if (proc_subset(0,0,0,0,is)) then

        ! fill zfshear_ga with the second derivative of phi_ga
        ! factor 0.5 due to normalisation with thermal velocity
        local_x_line_ga = 0.0
        do ipar = 1, ns
        do ix = 1, nx
          do j = 1, nmu
            ! put second radial derivatve of phi_ga in dum
            call d2fielddx2(PHI_GA_FIELD,iyzero,ipar,j,lsp(is),dum)
            do k = 1, nvpar
              ! the veloctiy space volumen element
              d3v = bn(ix,ipar)*intvp(ipar,j,k,lsp(is))*intmu(j)
              
              local_x_line_ga(ix) = local_x_line_ga(ix) + 0.5 &
                                     & * real(real(dum(ix))) * ints(ipar)  &
                                     & * d3v * fmaxwl(ix,ipar,j,k,lsp(is))
            end do
          end do
        end do 
        end do
      
        ! complete the flux surface average and the velocity space integral
        call mpireduce_sum_inplace(local_x_line_ga, shape(local_x_line_ga), &
           & COMM_SP_EQ_X_EQ)
      endif

      if(proc_subset(0,1,1,1,is) .or. root_processor) then
        !call gather_array(zfshear_ga_G, n_x_grid, number_of_species, &
        ! local_x_line_ga, nx, nsp, COMM_X_NE, COMM_SP_NE, .false.)
        call gather_array(zfshear_ga_G, local_x_line_ga, &
             mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
      endif

      
      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zfshear_ga(is), zfshear_ga_G, xy_fmt, ascii_fmt)
      endif


      if (proc_subset(0,0,0,0,is)) then
      
        ! fill zfflow_ga with the first derivative of phi_ga
        ! factor 0.5 due to normalisation with thermal velocity
        local_x_line_ga = 0.0
        do ipar = 1, ns
        do ix = 1, nx
          do j = 1, nmu
            ! put first radial derivatve of phi_ga in dum
            call dfielddx(PHI_GA_FIELD,iyzero,ipar,j,lsp(is),dum)
            do k = 1, nvpar
              ! the veloctiy space volumen element
              d3v = bn(ix,ipar)*intvp(ipar,j,k,lsp(is))*intmu(j)
              
              local_x_line_ga(ix) = local_x_line_ga(ix) + 0.5 &
                                     & * real(dum(ix)) * ints(ipar)  &
                                     & * d3v * fmaxwl(ix,ipar,j,k,lsp(is))
            end do
          end do
        end do 
        end do
 
        ! complete the flux surface average and the velocity space integral
        call mpireduce_sum_inplace(local_x_line_ga, shape(local_x_line_ga), &
         & COMM_SP_EQ_X_EQ)
    
      endif
    
      if(proc_subset(0,1,1,1,is) .or. root_processor) then
        !call gather_array(zfshear_ga_G, n_x_grid, number_of_species, &
        ! local_x_line_ga, nx, nsp, COMM_X_NE, COMM_SP_NE, .false.)
        call gather_array(zfflow_ga_G, local_x_line_ga, &
             mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
      endif
      
      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zfflow_ga(is), zfflow_ga_G, xy_fmt, ascii_fmt)
      endif
      
      
      if (proc_subset(0,0,0,0,is)) then
      
        ! fill zfphi_ga with phi_ga (without derivative)
        local_x_line_ga = 0.0
        do ipar = 1, ns
        do ix = 1, nx
          do j = 1, nmu
            do k = 1, nvpar
              d3v = bn(ix,ipar)*intvp(ipar,j,k,lsp(is))*intmu(j)
              local_x_line_ga(ix) = local_x_line_ga(ix) + &
                 & real(get_averaged_phi(iyzero,ix,ipar,j,lsp(is),fdisi)) &
                 & * ints(ipar) * d3v * fmaxwl(ix,ipar,j,k,lsp(is))
            end do
          end do
        end do
        end do
   
        ! complete the flux surface average and the velocity space integral
        call mpireduce_sum_inplace(local_x_line_ga, shape(local_x_line_ga), &
           & COMM_SP_EQ_X_EQ)
      endif
     
      if(proc_subset(0,1,1,1,is) .or. root_processor) then
        !call gather_array(zfshear_ga_G, n_x_grid, number_of_species, &
        ! local_x_line_ga, nx, nsp, COMM_X_NE, COMM_SP_NE, .false.)
        call gather_array(zfphi_ga_G, local_x_line_ga, &
             mpi_dtype_radial_line, COMM_X_NE, .true., is <= nsp)
      endif
      
      ! write on file 
      if (root_processor) then 
        call append_chunk(i_zfphi_ga(is), zfphi_ga_G, xy_fmt, ascii_fmt)
      endif
    
    end do loop_global_species
  
  
  
    !non-gyro-averaged phi
    if(proc_subset(0,0,1,1,1)) then
      
      !fill zfshear with the second derivative of phi times 
      ! factor 0.5 due to normalisation with thermal velocity
      local_x_line = 0.0
      do ipar = 1, ns
        call d2fielddx2(PHI_FIELD,iyzero,ipar,j,is,dum)
        do ix = 1, nx
          local_x_line(ix) = local_x_line(ix) + 0.5 &
                                & * real(dum(ix)) * ints(ipar)
        end do
      end do
      
      ! complete the flux surface avg
      call mpireduce_sum_inplace(local_x_line, shape(local_x_line), &
           & COMM_S_NE) 
      if(proc_subset(0,1,1,1,1)) then
        call gather_array(zfshear_G, n_x_grid, &
             & local_x_line, nx, COMM_X_NE, .false.)
      endif
        
        
      ! fill zfflow with the first derivative of phi
      local_x_line = 0.0
      do ipar = 1, ns
        call dfielddx(PHI_FIELD,iyzero,ipar,j,is,dum)
        do ix = 1, nx
          local_x_line(ix) = local_x_line(ix) + 0.5 &
                                & * real(dum(ix)) * ints(ipar)
        end do
      end do
        
      ! complete the flux surface avg
      call mpireduce_sum_inplace(local_x_line, shape(local_x_line), &
           & COMM_S_NE)
      if(proc_subset(0,1,1,1,1)) then
        call gather_array(zfflow_G, n_x_grid, &
             & local_x_line, nx, COMM_X_NE, .false.)
      endif
        
        
        
      ! fill zfphi with phi (without derivative)
      local_x_line = 0.0
      do ix = 1, nx
        do ipar = 1, ns
            local_x_line(ix) = local_x_line(ix) + &
               & real(phi(iyzero,ix,ipar)) * ints(ipar) 
        end do
      end do
     
      ! complete the flux surface avg
      call mpireduce_sum_inplace(local_x_line, shape(local_x_line), &
           & COMM_S_NE)
      if(proc_subset(0,1,1,1,1)) then
        call gather_array(zfphi_G, n_x_grid, &
             & local_x_line, nx, COMM_X_NE, .false.)
      end if
      
    end if    
    
  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()
    use io, only : append_chunk, xy_fmt, ascii_fmt
    use mpiinterface, only : root_processor

    if (.not.lcalc_zfshear) return

    call calc_largestep()

    if(root_processor) then
      !call append_chunk(i_zfshear_ga, zfshear_ga_G, xy_fmt, ascii_fmt)
      !call append_chunk(i_zfflow_ga, zfflow_ga_G, xy_fmt, ascii_fmt)
      !all append_chunk(i_zfphi_ga, zfphi_ga_G, xy_fmt, ascii_fmt)

      call append_chunk(i_zfshear, zfshear_G, xy_fmt, ascii_fmt)
      call append_chunk(i_zfflow, zfflow_G, xy_fmt, ascii_fmt)
      call append_chunk(i_zfphi, zfphi_G, xy_fmt, ascii_fmt)
      
    end if

    ! testing: compute the derivative of simple polynomials
    if(test_x_derivatives) then
      call output_test_x_derivatives()
    end if
  end subroutine output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_test_x_derivatives()
    use io, only : output_array, xy_fmt, ascii_fmt
    use grid, only : nx, n_x_grid, proc_subset
    use mpiinterface, only : gather_array, root_processor
    use mpicomms, only : COMM_X_NE
    use diagnos_generic, only : d2fielddx2
    use geom, only : dxgr
    
    if(proc_subset(0,1,1,1,1)) then
      ! uncomment this or that
      !call dfielddx(9999,1,1,1,1,local_x_line_ga)
      call d2fielddx2(9999,1,1,1,1,dum)
      
      call gather_array(zfflow_ga_G(:), n_x_grid, &
         & real(dum(:)), nx, COMM_X_NE, .false.)
      if(root_processor) then
        ! for testing radial derivatives:
        call output_array('x_derivs_test', '/diagnostic/diagnos_zfshear', &
           & zfflow_ga_G(:), 'F', xy_fmt, ascii_fmt)
        write(*,*) 'dxgr:', dxgr
      end if
    end if
  end subroutine output_test_x_derivatives

end module diagnos_zfshear
