!------------------------------------------------------------------------------
!> Computes and outputs fluxes in velocity space or full 5d.  
!> Should be consolidated with module fluxes, calculation presently duplicated.
!------------------------------------------------------------------------------
module diagnos_fluxes_vspace

  implicit none

  private

#define OLD_FLUX_OUTPUT

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  public :: calc_fluxes_full_detail
  public :: pflux_det, eflux_det, vflux_det
  public :: allocate_flux_det

  !> switch for output of the particle, momentum and heat fluxes on
  !> the 6D grid
  logical, save, public :: lfluxes_detail

  !> switches for output of the particle, momentum and heat fluxes on
  !> velspace, due to any field
  logical, save, public :: lfluxes_vspace
  !> switches for output of the particle, momentum and heat fluxes on
  !> velspace, due to any field
  logical, save, public :: lfluxes_vspace_phi
  logical, save, public :: lfluxes_vspace_em
  logical, save, public :: lfluxes_vspace_bpar

  logical, parameter :: lfluxes_detail_estep = .false.

  !> the particle flux  pflux(nmod,nx,number_of_species,ns,nmu,nvpar)
  !> due to the ExB motion
  real, target, save, allocatable, dimension(:,:,:,:,:,:) :: pflux_det

  !> the energy flux  eflux(nmod,nx,number_of_species,ns,nmu,nvpar)
  !> due to the ExB motion
  real, target, save, allocatable, dimension(:,:,:,:,:,:) :: eflux_det

  !> the toroidal angular momentum flux
  !> vflux(nmod,nx,number_of_species,ns,nmu,nvpar)
  !> due to the ExB motion
  real, target, save, allocatable, dimension(:,:,:,:,:,:) :: vflux_det

  !> the fluxes  [p/e/v]flux(number_of_species,nmu,nvpar) due
  !> to the ExB motion 
  real, save, allocatable, dimension(:,:,:) :: pflux_vspace, eflux_vspace
  real, save, allocatable, dimension(:,:,:) :: vflux_vspace

  ! the derivation of the fields with respect to s
  complex, allocatable :: dphi_gads(:), dapar_gads(:), dbpar_gads(:)

  integer, save :: mpi_dtype_vspace

  !> the range of MPI tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl


contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    lfluxes_detail = .false.
    lfluxes_vspace = .false.
    lfluxes_vspace_phi = .false.
    lfluxes_vspace_em = .false.
    lfluxes_vspace_bpar = .false. 
  end subroutine set_default_nml_values
  
  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(lfluxes_detail,1)
    call mpibcast(lfluxes_vspace,1)
    call mpibcast(lfluxes_vspace_phi,1)
    call mpibcast(lfluxes_vspace_em,1)
    call mpibcast(lfluxes_vspace_bpar,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine check()

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : MPIREAL_X, register_tag_range
    use mpidatatypes, only : create_subarray_datatype
    use global, only : id_vpar, id_mu, id_s, id_x, id_mod, id_sp
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD, DISTRIBUTION
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS
    use control, only : nlphi, nlapar, nlbpar
    use rho_par_switch, only : lflux_rhostar
    use grid, only : number_of_species
    logical, intent(inout) :: requirements(:,:)

    if (.not. lfluxes_detail) return

    requirements(DISTRIBUTION,LOCAL_DATA) = .true.

    if(lflux_rhostar) then
      if(nlphi) then
        requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
        requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
      if(nlapar) then
        requirements(APAR_GA_FIELD,LOCAL_DATA) = .true.
        requirements(APAR_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
      if(nlbpar) then
        requirements(BPAR_GA_FIELD,LOCAL_DATA) = .true.
        requirements(BPAR_GA_FIELD,S_GHOSTCELLS) = .true.
      end if
    end if    

    call create_subarray_datatype(MPIREAL_X,mpi_dtype_vspace, & 
       & id_vpar,id_mu,id_s,id_x,id_mod,id_sp)

    call register_tag_range(number_of_species * 3, &
         & tag_range_start, tag_range_end_inkl)

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use grid, only : nvpar, nmu, nsp
#ifndef OLD_FLUX_OUTPUT
    use mpiinterface, only : root_processor
    use grid, only : n_vpar_grid,n_mu_grid,n_s_grid,n_x_grid,nmod
    use grid, only : number_of_species
#endif
    
    use general, only : gkw_abort
    integer :: ierr

    call allocate_flux_det
    
    if (.not. (lfluxes_detail .or. &
       & lfluxes_vspace .or. &
       & lfluxes_vspace_phi .or. lfluxes_vspace_em .or. &
       & lfluxes_vspace_bpar)) return
    
#ifndef OLD_FLUX_OUTPUT
    if(lfluxes_detail_estep) then
      if(root_processor) then
        allocate(flux_det_global(n_vpar_grid,n_mu_grid,n_s_grid,n_x_grid,nmod,&
           & number_of_species),stat=ierr)
      else
        allocate(flux_det_global(1,1,1,1,1,1),stat=ierr)
      end if
      if (ierr /= 0) call gkw_abort('diagnostic :: flux_det_global')
    end if
#endif
    
    if(lfluxes_vspace .or. &
       & lfluxes_vspace_phi .or. lfluxes_vspace_em .or. &
       & lfluxes_vspace_bpar) then
      allocate(pflux_vspace(nvpar,nmu,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_vspace')
      allocate(eflux_vspace(nvpar,nmu,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_vspace')
      allocate(vflux_vspace(nvpar,nmu,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_vspace')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_flux_det
    use rho_par_switch, only : lflux_rhostar
    use grid, only : nmod,nx,nsp,ns,nmu,nvpar
    use general, only : gkw_abort
    use control, only : nlphi, nlapar, nlbpar
    
    integer :: ierr
    if(.not.allocated(pflux_det)) then
      allocate(pflux_det(nvpar,nmu,ns,nx,nmod,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: pflux_det')
      allocate(eflux_det(nvpar,nmu,ns,nx,nmod,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: eflux_det')
      allocate(vflux_det(nvpar,nmu,ns,nx,nmod,nsp),stat=ierr)
      if (ierr /= 0) call gkw_abort('diagnostic :: vflux_det')
      
    end if
    ! rho* correction of fluxes requires s-derivative of the fields
    if (lflux_rhostar) then
      if (nlphi .and..not. allocated(dphi_gads)) allocate(dphi_gads(ns))
      if (nlapar .and..not. allocated(dapar_gads)) allocate(dapar_gads(ns))
      if (nlbpar .and..not. allocated(dbpar_gads)) allocate(dbpar_gads(ns))
    endif

  end subroutine allocate_flux_det

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use mpidatatypes, only : free_datatype
    use rho_par_switch, only : lflux_rhostar
    
    if (.not. lfluxes_detail) return

    if (allocated(pflux_det)) deallocate(pflux_det)
    if (allocated(eflux_det)) deallocate(eflux_det)
    if (allocated(vflux_det)) deallocate(vflux_det)

    if (allocated(pflux_vspace)) deallocate(pflux_vspace)
    if (allocated(eflux_vspace)) deallocate(eflux_vspace)
    if (allocated(vflux_vspace)) deallocate(vflux_vspace)
    
#ifndef OLD_FLUX_OUTPUT
    if (allocated(flux_det_global)) deallocate(flux_det_global)
#endif

    if (lflux_rhostar)  then
      if(allocated(dphi_gads)) deallocate(dphi_gads)
      if(allocated(dapar_gads)) deallocate(dapar_gads)
      if(allocated(dbpar_gads)) deallocate(dbpar_gads)
    endif

    call free_datatype(mpi_dtype_vspace)

  end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output()
    use global, only : EVERY_FIELD
    
    if (lfluxes_detail) then
      call calc_fluxes_full_detail(EVERY_FIELD, .false.)
      call output_fluxes_full_detail(EVERY_FIELD)
    end if
  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD, EVERY_FIELD

    ! detailed velocity space fluxes
    if (lfluxes_vspace) then
      call calc_fluxes_full_detail(EVERY_FIELD, .false.)
      if(lfluxes_detail_estep) call output_fluxes_full_detail(EVERY_FIELD)
      call calc_fluxes_vspace
      call output_fluxes_vspace(EVERY_FIELD)
    end if

    if (lfluxes_vspace_phi) then
      call calc_fluxes_full_detail(PHI_FIELD, .false.)
      if(lfluxes_detail_estep) call output_fluxes_full_detail(PHI_FIELD)
      call calc_fluxes_vspace
      call output_fluxes_vspace(PHI_FIELD)
    end if
    if (lfluxes_vspace_em) then
      call calc_fluxes_full_detail(APAR_FIELD, .false.)
      if(lfluxes_detail_estep) call output_fluxes_full_detail(APAR_FIELD)
      call calc_fluxes_vspace
      call output_fluxes_vspace(APAR_FIELD)
    end if
    if (lfluxes_vspace_bpar) then
      call calc_fluxes_full_detail(BPAR_FIELD, .false.)
      if(lfluxes_detail_estep) call output_fluxes_full_detail(BPAR_FIELD)
      call calc_fluxes_vspace
      call output_fluxes_vspace(BPAR_FIELD)
    end if

  end subroutine output

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine calc_fluxes_vspace
    use mpiinterface, only : mpireduce_sum_inplace
    use grid, only : nsp, nmod, nx, ns, nmu, nvpar
    use mpicomms, only : COMM_S_NE_X_NE
    integer :: imod, ix, i, j, k, is
    pflux_vspace=0.
    eflux_vspace=0.
    vflux_vspace=0.

    !  Integrate
    do is = 1, nsp
      do imod = 1, nmod
        do ix = 1, nx
          do i = 1, ns

            ! velocity space
            do j = 1, nmu
              do k = 1, nvpar

                pflux_vspace(k,j,is) = pflux_vspace(k,j,is) + pflux_det(k,j,i,ix,imod,is)
                eflux_vspace(k,j,is) = eflux_vspace(k,j,is) + eflux_det(k,j,i,ix,imod,is)
                vflux_vspace(k,j,is) = vflux_vspace(k,j,is) + vflux_det(k,j,i,ix,imod,is)
              end do
            end do
          end do
        end do
      end do
    end do
    !complete the reductions over the s/x grid
    call mpireduce_sum_inplace(pflux_vspace,shape(pflux_vspace), COMM_S_NE_X_NE)
    call mpireduce_sum_inplace(eflux_vspace,shape(eflux_vspace), COMM_S_NE_X_NE)
    call mpireduce_sum_inplace(vflux_vspace,shape(vflux_vspace), COMM_S_NE_X_NE)

  end subroutine calc_fluxes_vspace
  
  !---------------------------------------------------------------------------
  !> This routine calculates and outputs the "decomposed" fluxes of particles,
  !> energy and parallel momentum, i.e. the integrals over s, mu and vpar are
  !> not performed.
  !>
  !> If the integral is performed, the total fluxes obtained are
  !>   flux = sum flux_det
  !>
  !> Note: the neoclassical part is not included.
  !---------------------------------------------------------------------------
  subroutine calc_fluxes_full_detail(field, magnitude)
    use grid,           only : nx, ns, nmu, nvpar, nsp, nmod
    use dist,           only : fdisi 
    use geom,           only : bn, efun, Rfun, signB, bt_frac,ints
    use mode,           only : krho, kxrh
    use components,     only : tmp, vthrat, signz, rhostar
    use velocitygrid,   only : intmu, intvp, mugr, vpgr
    use matdat,         only : get_f_from_g
    use fields,         only : get_averaged_phi, get_averaged_apar 
    use fields,         only : get_averaged_bpar
    use constants,      only : ci1
    use rho_par_switch, only : lflux_rhostar
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD
    use global, only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use diagnos_generic, only : dfieldds, parseval_correction
    use global,          only : gkw_a_equal_b_accuracy
    use control, only : nlphi, nlapar, nlbpar

    !> According to this argument, contributions caused by all or only a
    !> specific one of the field fluctuations are considered,
    integer, intent(in) :: field
    logical, intent(in) :: magnitude
    integer :: imod, ix, i, j, k, is

    ! Dummy variables 
    complex :: dum, dum1, dum2, fdis, phi_ga, apar_ga, bpar_ga
    complex :: dumes_rs, dumem_rs, dumbpar_rs

    ! switches for the contributions of the respective fields
    real    :: switch_phi, switch_apa, switch_bpa

    real :: d2X,d3X,d3v

    !switch off some contributions to the flux if required
    switch_phi=1.0
    switch_apa=1.0
    switch_bpa=1.0
    select case(field)
    case(PHI_FIELD)
      switch_apa=0.0
      switch_bpa=0.0
      if(.not.nlphi) then
        pflux_det = 0.0
        vflux_det = 0.0
        eflux_det = 0.0
        return
      end if
    case(APAR_FIELD)
      switch_phi=0.0
      switch_bpa=0.0
      if(.not.nlapar) then
        pflux_det = 0.0
        vflux_det = 0.0
        eflux_det = 0.0
        return
      end if
    case(BPAR_FIELD)
      switch_phi=0.0
      switch_apa=0.0
      if(.not.nlbpar) then
        pflux_det = 0.0
        vflux_det = 0.0
        eflux_det = 0.0
        return
      end if
    end select

    !  Calculate the fluxes
    d2X = 1.0
    do is = 1, nsp ; do imod = 1, nmod ; do ix = 1, nx ; do i = 1, ns
      d3X = ints(i) * d2X
      ! velocity space
      do j = 1, nmu

        ! rho* correction requires derivative of the fields 
        if (lflux_rhostar ) then
          select case(field)
          case(PHI_FIELD)
            call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dphi_gads)
          case(APAR_FIELD)
            call dfieldds(APAR_GA_FIELD,imod,ix,j,is,dapar_gads)
          case(BPAR_FIELD)  
            call dfieldds(BPAR_GA_FIELD,imod,ix,j,is,dbpar_gads)
          end select
        endif
        
        do k = 1, nvpar
          ! bn is the magnetic field modulus. The 2*pi which is
          ! furthermore contained in the velocity-space Jacobian is
          ! defined into intmu.
          d3v = intmu(j)*bn(ix,i) * intvp(i,j,k,is)

          ! the gyro-averaged fields 
          phi_ga  = switch_phi*get_averaged_phi(imod,ix,i,j,is,fdisi)
          apar_ga = switch_apa*get_averaged_apar(imod,ix,i,j,is,fdisi)
          bpar_ga = switch_bpa*get_averaged_bpar(imod,ix,i,j,is,fdisi)

          ! fdis is the distribution without A_par contribution
          fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi)

          ! If assuming 90 degree phase angle in fluxes,
          ! take the magnitudes of the complex numbers  
          if (magnitude) then
            phi_ga  = abs(phi_ga)
            apar_ga = abs(apar_ga) 
            bpar_ga = abs(bpar_ga)
            fdis    = -ci1*abs(fdis)
          end if

          ! in the implicit scheme fdis can be NaN for intvp = 0
          if (gkw_a_equal_b_accuracy(intvp(i,j,k,is), 0.0)) fdis = 0.

          dum = parseval_correction(imod)*(efun(ix,i,1,1)*kxrh(ix) + &
             & efun(ix,i,2,1)*krho(imod))*fdis

          dum1 = dum*(conjg(phi_ga)&
             & -2.*vthrat(is)*vpgr(i,j,k,is)*conjg(apar_ga)&
             & +2.*mugr(j)*tmp(ix,is)*conjg(bpar_ga)/signz(is))

          dum2 = dum*(conjg(phi_ga)*bn(ix,i)&
             & -2.*vthrat(is)*vpgr(i,j,k,is)*bn(ix,i)*conjg(apar_ga)&
             & +2.*mugr(j)*tmp(ix,is)*conjg(bpar_ga)/signz(is)*bn(ix,i))
          ! or simply
          !dum2 = dum1 * bn(ix,i)

          if (lflux_rhostar) then
            select case(field)
            case(PHI_FIELD)
              dumes_rs = parseval_correction(imod)&
                 & *efun(ix,i,3,1) * fdis &
                 * conjg(dphi_gads(i))
              !the fluxes requires the real part of the rho*
              !correction. below the imaginary part of dumes[12]
              !is used. therefore I multiply with ci1
              dum1 = dum1 + dumes_rs * rhostar * ci1
              dum2 = dum2 + dumes_rs*bn(ix,i)*rhostar * ci1
            case(APAR_FIELD)
              dumem_rs = -parseval_correction(imod)&
                 & *efun(ix,i,1,3)*bn(ix,i) * fdis * &
                 & vthrat(is)*vpgr(i,j,k,is) * conjg(dapar_gads(i))
                 
              dum1 = dum1 + rhostar * dumem_rs * ci1
              dum2 = dum2 + dumem_rs * bn(ix,i) * rhostar * ci1
            case(BPAR_FIELD)
              ! there was no bn in dumbpar_rs, in contrast to the
              ! other two _rs variables
              dumbpar_rs = parseval_correction(imod)&
                 & *efun(ix,i,1,3) * fdis  *                 &
                 & mugr(j)* tmp(ix,is)*conjg(dbpar_gads(i))&
                 & /(signz(is)*bn(ix,i))
              dum1 = dum1 + rhostar * dumbpar_rs * ci1
              dum2 = dum2 + dumbpar_rs*bn(ix,i)*rhostar * ci1
            end select
          end if

          pflux_det(k,j,i,ix,imod,is) = d3X*d3v*aimag(dum1)
          eflux_det(k,j,i,ix,imod,is) = d3X*d3v*(vpgr(i,j,k,is)**2*aimag(dum1) +&
             & 2.*mugr(j)*aimag(dum2))
          vflux_det(k,j,i,ix,imod,is) = d3X*d3v*(aimag(dum1)*vpgr(i,j,k,is)*Rfun(ix,i)* &
             & bt_frac(ix,i)*signB)

        end do
      end do

    end do; end do ; end do ; end do 

  end subroutine calc_fluxes_full_detail

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_fluxes_vspace(field)
    use io, only : binary_fmt
    use global, only : int2char_zeros
    use control, only : io_legacy
    use diagnos_generic, only : velocity_slice_output
    use grid, only : proc_subset, number_of_species, lsp, nsp
    use mpiinterface, only : root_processor
    integer, intent(in) :: field

    character (len=64) :: filename
    ! The actual species index is in isglb
    integer :: isglb, isp

    ! then gather and output in velocity space
    do isglb = 1, number_of_species
      
      if (proc_subset(1,1,0,0,isglb).or.root_processor) then

        filename="pfluxes_vspace"//trim(int2char_zeros(isglb,2))
        filename = trim(filename)//get_fieldname(field)

        if(io_legacy) then
          filename=trim(filename)//'_bin'
        end if
        ! filename=trim(filename)//"sp"//trim(file_count_suffix)
        ! To change these outputs to ASCII switch the line above 
        ! for the comment and remove. .true. in the calls below

        if(root_processor .and. isglb > nsp) then
          isp = 1
        else
          isp = lsp(isglb)
        end if

        call velocity_slice_output('diagnostic/diagnos_fluxes_vspace', &
           & pflux_vspace(:,:,isp),filename,1,1,isglb, binary_fmt, &
           & tag_range_start + (isglb-1)*3 + 0)

        filename(1:1)='e'
        call velocity_slice_output('diagnostic/diagnos_fluxes_vspace', &
           & eflux_vspace(:,:,isp),filename,1,1,isglb, binary_fmt, &
           & tag_range_start + (isglb-1)*3 + 1)

        filename(1:1)='v'
        call velocity_slice_output('diagnostic/diagnos_fluxes_vspace', &
           & vflux_vspace(:,:,isp),filename,1,1,isglb, binary_fmt, &
           & tag_range_start + (isglb-1)*3 + 2)
      end if
    end do

  end subroutine output_fluxes_vspace

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  pure function get_fieldname(field) result(fieldname)
    use global, only : PHI_FIELD, APAR_FIELD, BPAR_FIELD, EVERY_FIELD
    integer, intent(in) :: field
    character (len=4) :: fieldname
    select case(field)
    case(PHI_FIELD)
      fieldname = '_phi'
    case(APAR_FIELD)
      fieldname = '_apa'
    case(BPAR_FIELD)
      fieldname = '_bpa'
    case(EVERY_FIELD)
      fieldname = ''
    end select
  end function get_fieldname

  !---------------------------------------------------------------------------
  !>
  !---------------------------------------------------------------------------
  subroutine output_fluxes_full_detail(field)
#ifdef OLD_FLUX_OUTPUT
    use io_binary, only : mpifile_open, mpifile_close
    use mpiinterface, only : MPI_OFFSET_KIND, MPI_INFO_NULL
    use mpiinterface, only : MPIREAL_X, statusmpi
#else
    use control, only : time
    use io, only : attach_metadata, binary_fmt, xy_fmt, mpi_output_array
    use io, only : comments_key, description_key, phys_unit_key
    use io, only : not_avail
    use grid, only : proc_subset
#endif
    use mpicomms, only : COMM_CART
    integer, intent(in) :: field
#if defined(mpi2)
#ifdef OLD_FLUX_OUTPUT
    integer :: file_unit, ierr
    integer(kind=MPI_OFFSET_KIND), parameter :: view_disp = 0
#endif
    character (len=64) :: filename
    
    filename = "fluxes_det"//trim(get_fieldname(field))//".dat"

#ifdef OLD_FLUX_OUTPUT
    ! If one uses that older code block, by defining OLD_FLUX_OUTPUT,
    ! one will have binary MPI-IO for any io_format.

    call mpifile_open(trim(filename),file_unit,COMM_CART,.true.)
    call mpi_file_set_view(file_unit,view_disp,MPIREAL_X,mpi_dtype_vspace, &
       & 'native', MPI_INFO_NULL,ierr)
    call mpi_file_write_all(file_unit,real(pflux_det),size(pflux_det), &
       & MPIREAL_X,statusmpi,ierr)
    call mpi_file_write_all(file_unit,real(eflux_det),size(eflux_det), &
       & MPIREAL_X,statusmpi,ierr)
    call mpi_file_write_all(file_unit,real(vflux_det),size(vflux_det), &
       & MPIREAL_X,statusmpi,ierr)
    call mpifile_close(file_unit)
#else

    ! If one uses that newer code block, by undefining OLD_FLUX_OUTPUT,
    ! one will have binary MPI-IO only for io_format='binary' or 'mixed'.

    ! If one day in the far future we are able to use PHDF5, we can
    ! use that code block permanently. For now it sketches at least
    ! the possibility to have them in hdf5, for not too large runs:
    ! Output is serial and will involve a 6d gather communication for
    ! ascii and hdf5 format, which is a bad idea for large runs.
    filename = 'pflux_detail'
    call mpi_output_array(filename, 'diagnostic/diagnos_fluxes_vspace', &
       & real(pflux_det), mpi_dtype_vspace, &
       & flux_det_global, &
       & COMM_CART, xy_fmt, binary_fmt, proc_subset(0,0,0,0,0), &
       & .true.)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', phys_unit_key, not_avail, &
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', description_key, &
       & '',&
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', comments_key, &
       & 'ints and bn*intmu*intvp are already included.', binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', 'time', time, binary_fmt)

    filename = 'eflux_detail'
    call mpi_output_array(filename, 'diagnostic/diagnos_fluxes_vspace', &
       & real(eflux_det), mpi_dtype_vspace, &
       & flux_det_global, &
       & COMM_CART, xy_fmt, binary_fmt, proc_subset(0,0,0,0,0), &
       & .true.)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', phys_unit_key, not_avail, &
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', description_key, &
       & '',&
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', comments_key, &
       & 'ints and bn*intmu*intvp are already included.', binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', 'time', time, binary_fmt)

    filename = 'vflux_detail'
    call mpi_output_array(filename, 'diagnostic/diagnos_fluxes_vspace', &
       & real(vflux_det), mpi_dtype_vspace, &
       & flux_det_global, &
       & COMM_CART, xy_fmt, binary_fmt, proc_subset(0,0,0,0,0), &
       & .true.)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', phys_unit_key, not_avail, &
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', description_key, &
       & '',&
       & binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', comments_key, &
       & 'ints and bn*intmu*intvp are already included.', binary_fmt)
    call attach_metadata(filename, &
       & 'diagnostic/diagnos_fluxes_vspace', 'time', time, binary_fmt)
#endif

#endif

  end subroutine output_fluxes_full_detail

end module diagnos_fluxes_vspace
