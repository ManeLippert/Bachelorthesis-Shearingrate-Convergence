!------------------------------------------------------------------------------
!>
!> 
!>
!------------------------------------------------------------------------------
module diagnos_jdote

  implicit none

  private
  
  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: calc_smallstep
  public :: finalize
  public :: output
  
  !> True if the flux surface averaged J.E diagnostic is switched on
  logical, save, public :: lcalc_jdote_fs

  !> True if the J.E diagnostic is switched on
  logical, save, public :: lcalc_jdote

  !> time between call of calc_smallstep() and output()
  real, save :: delta_time_jdote = 1.
  
  !> last apar field, for real calculation of the
  !> time derivative of apar.
  complex, save, allocatable, dimension(:,:,:,:,:) :: last_aparga

  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl
  integer, save :: tag_range_start2, tag_range_end_inkl2
  integer, parameter :: n_velocity_slices = 6
contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

    lcalc_jdote_fs = .false.
    lcalc_jdote = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_jdote,1)
    call mpibcast(lcalc_jdote_fs,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_abort, gkw_warn
    use mode, only : mode_box, lxn, lyn
    use diagnos_velspace, only : npointsvel, psi_velspace, zeta_velspace
    integer :: ipoint
    real :: psi, zeta
    
    if (lcalc_jdote_fs .and..not.mode_box) then
      call gkw_warn('the lcalc_jdote_fs diagnostic needs to do a Fourier trafo&
         &, but this is not possible for mode_box=F. Does this diagnostic &
         & really make sense then?')
      lcalc_jdote_fs = .false.
    end if

    if (lcalc_jdote .and..not.mode_box) then
      call gkw_warn('the lcalc_jdote diagnostic needs to do a Fourier trafo&
         &, but this is not possible for mode_box=F. Does this diagnostic &
         & really make sense then?')
      lcalc_jdote = .false.
    end if

    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return
    
    if (npointsvel <= 0) then
      call gkw_abort('lcalc_jdote=T is set, but npointsvel is 0. &
         & You have to specify positions to evaluate this diagnostic.')
    end if

    do ipoint=1,npointsvel
      psi=psi_velspace(ipoint)
      if (psi > lxn/2. .or. psi < -lxn/2.) then
        write(*,*) 'psi_velspace:', psi, 'lxn:', lxn
        call gkw_abort('bad psi_velspace point selection')
      end if

      if (lcalc_jdote) then
        zeta=zeta_velspace(ipoint)
        if (zeta > lyn/2. .or. zeta < -lyn/2.) then
          write(*,*) 'zeta_velspace:', zeta
          call gkw_abort('bad zeta_velspace point selection')
        end if
      end if
    end do

    ! Now this diagnostic should work also in global simulations
    !if (.not.flux_tube) then
    !  call gkw_abort('Cannot use diagnos_jdote in a global run, the correct &
    !     & dependence of T on x has never been implemented')
    !endif

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor, register_tag_range
    use global, only : DISTRIBUTION, PHI_GA_FIELD
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS, X_GHOSTCELLS
    use diagnos_velspace, only : npointsvel
    use grid, only : number_of_species
    logical, intent(inout) :: requirements(:,:)

    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return

    call register_tag_range(npointsvel*number_of_species*n_velocity_slices, &
       & tag_range_start, tag_range_end_inkl)
    call register_tag_range(npointsvel*number_of_species*n_velocity_slices, &
       & tag_range_start2, tag_range_end_inkl2)

    requirements(DISTRIBUTION,LOCAL_DATA) = .true.
    requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
    requirements(PHI_GA_FIELD,X_GHOSTCELLS) = .true.
    
    if(root_processor) then
      ! Open logical units here.

      ! For per-timestep quantities, open the logical unit only if
      ! lwrite_output1 is set.
    end if
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use general, only : gkw_abort
    use grid, only : nmod, nx, ns, nmu, nsp

    integer :: ierr

    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return
    
    ! allocate(local_array_name(nmod,number_of_species),stat=ierr)
    ! if (ierr /= 0) &
    !    & call gkw_abort('diagnos_jdote :: could not allocate array')
    allocate(last_aparga(nmod,nx,ns,nmu,nsp),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: last_aparga')

  end subroutine allocate_mem


  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()

    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return
    
    ! deallocate all arrays of this diagnostic
    !if(allocated(local_array_name)) deallocate(local_array_name)
    if(allocated(last_aparga)) deallocate(last_aparga)
  end subroutine finalize


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
    
    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return

  end subroutine calc_largestep

  subroutine calc_smallstep(i_smallstep)
    use control,          only : naverage, time, flux_tube, nlapar
    integer, intent(in) :: i_smallstep
    real, save :: last_time_pre_jdote = 0.
    real, save :: delta_time_jdote = 1.

    if ((.not.lcalc_jdote).and.(.not. lcalc_jdote_fs)) return
    if (.not.(flux_tube .and. nlapar)) return

    ! if (i_smallstep == n_smallstep - 1) then
    ! else if (i_smallstep == n_smallstep) then
    ! end if

    ! FIXME This calculation of a time derivative presumably does
    ! not work for naverage == 1

    if (i_smallstep == naverage - 1) then
      call calc_aparga(last_aparga)
      ! store the time this subroutine was last called in a local variable
      last_time_pre_jdote = time

    else if (i_smallstep == naverage) then
      ! save time interval between this step and the
      ! last call of calc_entr().
      ! this intervall is needed by several subroutines of this module
      delta_time_jdote = time - last_time_pre_jdote
    end if

  end subroutine calc_smallstep
  
  !----------------------------------------------------------------------------
  !>
  !----------------------------------------------------------------------------
  subroutine calc_aparga(aparga)
    use grid,           only : nmod, nx, ns, nmu, nsp
    use dist,           only : fdisi
    use fields,         only : get_averaged_apar
    complex, intent(out) :: aparga(nmod,nx,ns,nmu,nsp)
    integer :: i, j, is, imod, ix

    do imod = 1, nmod
      do ix = 1, nx
        do i = 1, ns
          do j = 1, nmu
            do is = 1, nsp
              aparga(imod,ix,i,j,is) = get_averaged_apar(imod,ix,i,j,is,fdisi)
            end do
          end do
        end do
      end do
    end do

  end subroutine calc_aparga


  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output(file_count)
    integer, intent(in) :: file_count

    if((.not. lcalc_jdote).and.(.not. lcalc_jdote_fs)) return

    if (lcalc_jdote) then
      call output_jdote_pointwise(file_count)
    end if
    if (lcalc_jdote_fs) then
      call output_jdote_fs_pointwise(file_count)
    end if

  end subroutine output

  !****************************************************************************
  !> Routine to output velocity space at selected (x,y) points.
  !> These points are selected using psi_velspace and zeta_velspace
  !> and npointsvel in the diagnostic namelist.  The values can go from 
  !> psi = -/+ LX/2 and zeta = -/+ LY/2 relative to the xphi and yphi grids
  !> where 0,0 is centred on the island center.
  !> The values are integrated along the parallel direction (can be changed)
  !> and the integrand (currently radial / poloidal drifts) can also be changed.
  !----------------------------------------------------------------------------
  subroutine output_jdote_pointwise(file_count)
    use mpiinterface,     only : mpiallreduce_sum, mpireduce_sum_inplace
    use mpicomms,         only : COMM_S_NE_X_NE
    use control,          only : flux_tube, nlapar, time
    use geom,             only : ints, bn, efun, ffun
    use mode,             only : lxn, lyn
    use velocitygrid,     only : intvp, intmu, vpgr
    use global,           only : int2char_zeros
    use grid,             only : ns, nmu, nvpar, n_x_grid, psil, psih
    use grid,             only : lsp, number_of_species, lrx, proc_subset
    use general,          only : gkw_abort 
    use dist,             only : get_apar
    use components,       only : signz, tearingmode
    use components,       only : rhostar, vthrat, signz
    use dist,             only : iphi
    use diagnos_generic,  only : velocity_slice_output, dphigads_xy_point_get
    use diagnos_generic,  only : dphigadx_xy_point_get, dphigadzeta_xy_point_get
    use diagnos_generic,  only : dapargads_xy_point_get, dapargadx_xy_point_get
    use diagnos_generic,  only : dapargadzeta_xy_point_get, phiga_xy_point_get
    use io,               only : binary_fmt
    use linear_terms,     only : drift
    use rotation,         only : coriolis, cf_drift
    use tearingmodes,     only : omega_rot
    use constants,        only : pi
    use diagnos_velspace, only : npointsvel, psi_velspace, zeta_velspace, fdis_xy_point_get

    integer, intent(in) :: file_count
    real :: psi, zeta, fd
    real :: phiga, dphigads, dphigadx, dphigadzeta
    real :: dapargads, dapargadx, dapargadzeta, dapargadt
    integer :: i,j,k,ispl,ispg, ipoint, ixg, ix
    character (len=64) :: luname

    real :: drift_x, drift_y, drift_z

    real, allocatable, save, dimension(:,:) :: &
      & jdote_para_elec, jdote_perp_elec, jdote,&
      & jdote_para_magnetic, jdote_em, ene_pot

    character (len=6), save :: file_count_suffix
    file_count_suffix=trim(int2char_zeros(file_count,6))

    if (npointsvel < 1) return

    if (.not.lcalc_jdote) return

    !allocate array to contain the local slice
    if (.not. allocated(jdote_para_elec)) allocate(jdote_para_elec(nvpar,nmu))
    if (.not. allocated(jdote_perp_elec)) allocate(jdote_perp_elec(nvpar,nmu))
    if (.not. allocated(jdote_para_magnetic)) allocate(jdote_para_magnetic(nvpar,nmu))
    if (.not. allocated(jdote)) allocate(jdote(nvpar,nmu))
    if (.not. allocated(jdote_em)) allocate(jdote_em(nvpar,nmu))
    if (.not. allocated(ene_pot)) allocate(ene_pot(nvpar,nmu))

    ! Loop over the output points you are interested in
    do ipoint=1,npointsvel

      !The coordinates of the points you are interested in have to
      !be introduced referring to the xphi-yphi grid AFTER the
      !reorganization of the matrix (see below)
      psi=psi_velspace(ipoint)
      zeta=zeta_velspace(ipoint)
      ! psi = x/rho in flux tube simulations
      ! psi = r/R in global simulations

      !Updating zeta for rotating islands
      if (tearingmode) then
        zeta=zeta+lyn/2+lyn*omega_rot*time/2/pi
      end if
      ! Zeta is put in such that zeta=0 corresponds to the O-point (at t=0)
      ! Psi is zero at the center of the box AFTER THE REORGANIZATION of the 
      ! matrices. Recall that, for magnetic islands, the q=m/n surface 
      ! (i.e. the symmetry axes of the island) is at the edge of the box, 
      ! rather than on the centre, and therefore the domain must be reorganized 
      ! (in post-processing) in order to find the O-point
      ! at the centre of the box. The value of psi_velspace and
      ! zeta_velspace to be introduced refer to the grid AFTER the reorganization
      ! Hence 0,0 refers to the centre of the island    

      !IMPORTANT: if you do not reorganize the matrix (i.e. no magnetic island), 
      !uncomment the following line:
      ! We use lxn since it is the width of the box normalised to rho_ref
      ! Also, psi is the radial position, normalised to rho_ref, in flux tube simulations
      ! In global simulations psi = r/R
      if (flux_tube) then
        psi=psi+lxn/2 ! This is to place psi=0 at the centre of the box
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if        

      ix = lrx(ixg)

      ! loop over global species
      do ispg = 1, number_of_species
        jdote_para_elec(:,:) = 0.0
        jdote_para_magnetic(:,:) = 0.0
        jdote_perp_elec(:,:) = 0.0
        jdote_em(:,:) = 0.0
        jdote(:,:) = 0.0
        ene_pot(:,:) = 0.0
        ! We need to determine the local species index
        ispl=lsp(ispg)
        !if that species is not on this processor skip to next species
        if (proc_subset(ixg,0,0,0,ispg)) then
          do j = 1, nmu; do k=1,nvpar
            do i=1,ns
             ! obtain the drift
              call drift(ix,i,j,k,ispl,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
              fd = fdis_xy_point_get(psi,zeta,i,j,k,ispl)
              dphigadx = dphigadx_xy_point_get(psi,zeta,i,j,ispl)
              dphigadzeta = dphigadzeta_xy_point_get(psi,zeta,i,j,ispl)
              dphigads = dphigads_xy_point_get(psi,zeta,i,j,ispl)
              ! This is the current \cdot E-field
              ! i.e. exchange of energy particles-waves
              jdote_perp_elec(k,j) = jdote_perp_elec(k,j) &
                & + signz(ispl) * (drift_x*dphigadx + drift_y*dphigadzeta) &
                & * fd*bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
              jdote_para_elec(k,j) = jdote_para_elec(k,j) &
                  & + signz(ispl) * vpgr(i,j,k,ispl)*vthrat(ispl)*ffun(ix,i)*dphigads &
                  & * fd*bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
              if (nlapar) then
                dapargadt = dapargadt_xy_point_get(psi,zeta,i,j,ispl)
                jdote_para_magnetic(k,j) = jdote_para_magnetic(k,j) &
                  & + 2.0*signz(ispl)*vpgr(i,j,k,ispl)*vthrat(ispl)*dapargadt &
                  & * fd*bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
                ! Calculation of v_\chi\cdot\nabla\phi
                dapargadx = dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                dapargadzeta = dapargadzeta_xy_point_get(psi,zeta,i,j,ispl)
                jdote_em(k,j) = jdote_em(k,j) + 2.0 * signz(ispl) * vpgr(i,j,k,ispl) * vthrat(ispl) * &
                  & efun(ix,i,1,2) * (dphigadx*dapargadzeta - dphigadzeta*dapargadx) * & ! psi,zeta components
                  & fd * bn(ix,i) * intvp(i,j,k,ispl) * intmu(j) * ints(i)
              end if
              ! FIRST TRY to extend to global simulations
              if (.not. flux_tube) then
                dphigads = dphigads_xy_point_get(psi,zeta,i,j,ispl)
                jdote_para_elec(k,j) = jdote_para_elec(k,j) &
                  & + signz(ispl) * rhostar*drift_z*dphigads &
                  & * fd*bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

                if (nlapar) then
                  dapargadt = dapargadt_xy_point_get(psi,zeta,i,j,ispl)
                  jdote_para_magnetic(k,j) = jdote_para_magnetic(k,j) &
                    & + signz(ispl) * rhostar*drift_z*dapargadt &
                    & * fd*bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
                  ! Calculation of v_\chi\cdot\nabla\phi
                  dapargads = dapargads_xy_point_get(psi,zeta,i,j,ispl)
                  dapargadx = dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                  dapargadzeta = dapargadzeta_xy_point_get(psi,zeta,i,j,ispl)
                  jdote_em(k,j) = jdote_em(k,j) + &
                    & signz(ispl) * rhostar * 2.0 * vpgr(i,j,k,ispl) * vthrat(ispl) * ( &
                    & efun(ix,i,1,3) * (dphigadx*dapargads - dphigads*dapargadx) + & ! psi,s components
                    & efun(ix,i,2,3) * (dphigadzeta*dapargads - dphigads*dapargadzeta))* & !zeta,s components
                    & fd * bn(ix,i) * intvp(i,j,k,ispl) * intmu(j) * ints(i)
                end if
              end if
              phiga = phiga_xy_point_get(psi,zeta,i,j,ispl)
              ene_pot(k,j) = ene_pot(k,j) + phiga * &
                & fd * bn(ix,i) * intvp(i,j,k,ispl) * intmu(j) * ints(i)
            end do  ! Loop over s (integration)
            jdote(k,j) = jdote_para_elec(k,j) + jdote_para_magnetic(k,j) + & 
              & jdote_perp_elec(k,j) + jdote_em(k,j)
          end do; end do ! Loop over vpar and mu
        else
          jdote(:,:) = 0.0
          jdote_para_elec(:,:) = 0.0
          jdote_para_magnetic(:,:) = 0.0
          jdote_perp_elec(:,:) = 0.0
          jdote_em(:,:) = 0.0
          ene_pot(:,:) = 0.0
        end if

        call mpireduce_sum_inplace(jdote_para_elec, shape(jdote_para_elec), & 
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_para_magnetic, shape(jdote_para_magnetic), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_perp_elec, shape(jdote_perp_elec), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_em, shape(jdote_em), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote, shape(jdote), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(ene_pot, shape(ene_pot), &
           & COMM_S_NE_X_NE)

        luname="jdote_para_elec_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_para_elec,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 0)

        luname="jdote_para_magnetic_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_para_magnetic,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 1)

        luname="jdote_perp_elec_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_perp_elec,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 2)

        luname="jdote_em_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_em,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 3)

        luname="jdote_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 4)

        luname="ene_pot_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & ene_pot,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 5)
      end do
    end do

  end subroutine output_jdote_pointwise

  !**********************************************************************
  !> Routine to output velocity space at selected x points.
  !> These points are selected using psi_velspace
  !> and npointsvel in the diagnostic namelist. The values can go from 
  !> psi = -/+ LX/2 relative to the xphi grid
  !> where 0,0 is centred on the island center.
  !> The values are integrated along the parallel and zeta directions
  !----------------------------------------------------------------------
  subroutine output_jdote_fs_pointwise(file_count)
    use mpiinterface,     only : mpiallreduce_sum, mpireduce_sum_inplace
    use mpicomms,         only : COMM_S_NE_X_NE
    use control,          only : flux_tube, nlapar, non_linear
    use geom,             only : ints, bn, efun, ffun
    use mode,             only : lxn, krho
    use velocitygrid,     only : intvp, intmu, vpgr
    use global,           only : int2char_zeros
    use grid,             only : ns, nmu, nvpar, n_x_grid, psil, psih, nmod
    use grid,             only : lsp, number_of_species, lrx, proc_subset, nx
    use general,          only : gkw_abort
    use dist,             only : get_apar, fdisi
    use components,       only : rhostar, vthrat, signz, dgrid
    use dist,             only : iphi, fdisi
    use diagnos_generic,  only : velocity_slice_output, dphigads_xy_point_get
    use diagnos_generic,  only : dfielddx, dfieldds
    use global,           only : PHI_GA_FIELD, APAR_GA_FIELD
    use io,               only : binary_fmt
    use linear_terms,     only : drift
    use rotation,         only : coriolis, cf_drift
    use constants,        only : pi
    use diagnos_velspace, only : npointsvel, psi_velspace
    use matdat,           only : get_f_from_g
    use constants,        only : ci1
    use fields,           only : get_averaged_phi, get_averaged_apar

    integer, intent(in) :: file_count
    real :: psi
    complex :: fd, dapargadt
    complex, allocatable :: dphigads(:), dphigadx(:)
    complex, allocatable :: dapargads_tmp(:), dapargadx_tmp(:)
    complex, allocatable :: dphigads_tmp(:), dphigadx_tmp(:)
    complex :: phiga, phiga_tmp, aparga_tmp
    complex :: vchi_gradphi
    integer :: i,j,k,ispl,ispg, ipoint, ixg, ix, imod, imod_1, imod_2
    character (len=64) :: luname

    real :: drift_x, drift_y, drift_z

    real, allocatable, save, dimension(:,:) :: &
      & jdote_para_elec, jdote_perp_elec, &
      & jdote, jdote_para_magnetic, &
      & jdote_em, jdote_perp_elec_rhostar, jdote_em_rhostar

    character (len=6), save :: file_count_suffix
    file_count_suffix=trim(int2char_zeros(file_count,6))

    if (npointsvel < 1) return

    if (.not.lcalc_jdote_fs) return

    !allocate array to contain the local slice
    if (.not. allocated(jdote_para_elec)) allocate(jdote_para_elec(nvpar,nmu))
    if (.not. allocated(jdote_perp_elec)) allocate(jdote_perp_elec(nvpar,nmu))
    if (.not. allocated(jdote_para_magnetic)) allocate(jdote_para_magnetic(nvpar,nmu))
    if (.not. allocated(jdote)) allocate(jdote(nvpar,nmu))
    if (.not. allocated(jdote_em)) allocate(jdote_em(nvpar,nmu))
    if (.not. allocated(jdote_perp_elec_rhostar)) allocate(jdote_perp_elec_rhostar(nvpar,nmu))
    if (.not. allocated(jdote_em_rhostar)) allocate(jdote_em_rhostar(nvpar,nmu))

    !allocate intermediate vectors for derivatives
    if(.not. allocated(dphigads)) allocate(dphigads(ns))
    if(.not. allocated(dphigadx)) allocate(dphigadx(nx))
    if(.not. allocated(dapargads_tmp)) allocate(dapargads_tmp(ns))
    if(.not. allocated(dapargadx_tmp)) allocate(dapargadx_tmp(nx))
    if(.not. allocated(dphigads_tmp)) allocate(dphigads_tmp(ns))
    if(.not. allocated(dphigadx_tmp)) allocate(dphigadx_tmp(nx))

    ! Loop over the output points you are interested in
    do ipoint=1,npointsvel

      !The coordinates of the points you are interested in have to
      !be introduced referring to the xphi-yphi grid AFTER the
      !reorganization of the matrix (see below)
      psi=psi_velspace(ipoint)
      ! psi = x/rho in flux tube simulations
      ! psi = r/R in global simulations

      !IMPORTANT: if you do not reorganize the matrix (i.e. no magnetic island), 
      !uncomment the following line:
      ! We use lxn since it is the width of the box normalised to rho_ref
      ! Also, psi is the radial position, normalised to rho_ref, in flux tube simulations
      ! In global simulations psi = r/R
      if (flux_tube) then
        psi=psi+lxn/2 ! This is to place psi=0 at the centre of the box
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if

      ix = lrx(ixg)

      ! loop over global species
      do ispg = 1, number_of_species
        jdote_para_elec(:,:) = 0.0
        jdote_para_magnetic(:,:) = 0.0
        jdote_perp_elec(:,:) = 0.0
        jdote_perp_elec_rhostar(:,:) = 0.0
        jdote_em(:,:) = 0.0
        jdote_em_rhostar(:,:) = 0.0
        jdote(:,:) = 0.0
        ! We need to determine the local species index
        ispl=lsp(ispg)

        !if that species is not on this processor skip to next species
        if (proc_subset(ixg,0,0,0,ispg)) then
          do j = 1, nmu; do k=1, nvpar
            do i=1, ns; do imod=1, nmod !Flux-surface average
              ! obtain the drift
              call drift(ix,i,j,k,ispl,drift_x,drift_y,drift_z,coriolis,cf_drift,.true.)
              fd = get_f_from_g(imod,ix,i,j,k,ispl,fdisi) ! fdis_xy_point_get(psi,zeta,i,j,k,ispl)
              call dfielddx(PHI_GA_FIELD,imod,i,j,ispl,dphigadx) !dphigadx_xy_point_get(psi,zeta,i,j,ispl)
              phiga = get_averaged_phi(imod,ix,i,j,ispl,fdisi) !dphigadzeta = dphigadzeta_xy_point_get(psi,zeta,i,j,ispl)
              call dfieldds(PHI_GA_FIELD,imod,ix,j,ispl,dphigads) !dphigads_xy_point_get(psi,zeta,i,j,ispl)
              ! This is the current \cdot E-field
              ! i.e. exchange of power particles-waves
              jdote_perp_elec(k,j) = jdote_perp_elec(k,j) &
                & + signz(ispl) &
                & * real((drift_x*conjg(dphigadx(ix)) + drift_y*conjg(ci1*krho(imod)*phiga)) &
                & * dgrid(ispl) * fd) &
                & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

              jdote_para_elec(k,j) = jdote_para_elec(k,j) &
                & + signz(ispl) &
                & * (vpgr(i,j,k,ispl)*vthrat(ispl)*ffun(ix,i)) &
                & * real(conjg(dphigads(i)) * dgrid(ispl) * fd) &
                & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

              if (nlapar) then
                dapargadt = dapargadt_x_point_get(psi,imod,i,j,ispl)
                jdote_para_magnetic(k,j) = jdote_para_magnetic(k,j) &
                  & + signz(ispl) &
                  & * (2.0*vpgr(i,j,k,ispl)*vthrat(ispl)) &
                  & * real(conjg(dapargadt) * dgrid(ispl) * fd) &
                  & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

                if (non_linear) then
                  ! Calculation of v_\chi\cdot\nabla\phi
                  vchi_gradphi = (0.0, 0.0)
                  do imod_1 = 1,nmod
                    imod_2 = imod - imod_1 + 1
                    if (imod_2.gt.0) then
                      call dfielddx(APAR_GA_FIELD,imod_1,i,j,ispl,dapargadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      !call dfieldds(APAR_GA_FIELD,imod_1,i,j,ispl,dapargads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      call dfielddx(PHI_GA_FIELD,imod_2,i,j,ispl,dphigadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      !call dfieldds(PHI_GA_FIELD,imod_2,i,j,ispl,dphigads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      phiga_tmp = get_averaged_phi(imod_2,ix,i,j,ispl,fdisi)
                      aparga_tmp = get_averaged_apar(imod_1,ix,i,j,ispl,fdisi)
                      vchi_gradphi = vchi_gradphi &
                        & + 2.0 * vpgr(i,j,k,ispl) * vthrat(ispl) * &
                        & efun(ix,i,1,2) * ((dphigadx_tmp(ix))*aparga_tmp*ci1*krho(imod_1) &
                        & - (ci1*krho(imod_2)*phiga_tmp)*dapargadx_tmp(ix)) ! psi, zeta components
                    else
                      imod_2 = imod_1 - imod + 1
                      call dfielddx(APAR_GA_FIELD,imod_1,i,j,ispl,dapargadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      !call dfieldds(APAR_GA_FIELD,imod_1,i,j,ispl,dapargads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      call dfielddx(PHI_GA_FIELD,imod_2,i,j,ispl,dphigadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      !call dfieldds(PHI_GA_FIELD,imod_2,i,j,ispl,dphigads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      phiga_tmp = get_averaged_phi(imod_2,ix,i,j,ispl,fdisi)
                      aparga_tmp = get_averaged_apar(imod_1,ix,i,j,ispl,fdisi)
                      vchi_gradphi = vchi_gradphi &
                        & + 2.0 * vpgr(i,j,k,ispl) * vthrat(ispl) * &
                        & efun(ix,i,1,2) * (conjg(dphigadx_tmp(ix))*aparga_tmp*ci1*krho(imod) &
                        & - conjg(ci1*krho(imod)*phiga_tmp)*dapargadx_tmp(ix)) ! psi, zeta components
                    end if
                    
                  end do
                  jdote_em(k,j) = jdote_em(k,j) &
                    & + signz(ispl) * real(conjg(vchi_gradphi) * dgrid(ispl) * fd) &
                    & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
                end if !non_linear
              end if !nlapar
              ! FIRST TRY to extend to global simulations
              if (.not. flux_tube) then
                jdote_perp_elec_rhostar(k,j) = jdote_perp_elec_rhostar(k,j) &
                  & + signz(ispl) &
                  & * (rhostar*drift_z) &
                  & * real(conjg(dphigads(i)) * dgrid(ispl) * fd) &
                  & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

                if (nlapar.and.non_linear) then
                  
                  ! Calculation of v_\chi\cdot\nabla\phi
                  vchi_gradphi = (0.0, 0.0)
                  do imod_1 = 1,nmod
                    imod_2 = imod - imod_1 + 1
                    if (imod_2.gt.0) then
                      call dfielddx(APAR_GA_FIELD,imod_1,i,j,ispl,dapargadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      call dfieldds(APAR_GA_FIELD,imod_1,i,j,ispl,dapargads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      call dfielddx(PHI_GA_FIELD,imod_2,i,j,ispl,dphigadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      call dfieldds(PHI_GA_FIELD,imod_2,i,j,ispl,dphigads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      phiga_tmp = get_averaged_phi(imod_2,ix,i,j,ispl,fdisi)
                      aparga_tmp = get_averaged_apar(imod_1,ix,i,j,ispl,fdisi)
                      vchi_gradphi = vchi_gradphi &
                        & + 2.0 * vpgr(i,j,k,ispl) * vthrat(ispl) * ( &                        
                        & efun(ix,i,1,3) * rhostar * ((dphigadx_tmp(ix))*dapargads_tmp(i) &
                          & - (dphigads_tmp(i))*dapargadx_tmp(ix)) + & ! psi, s components
                        & efun(ix,i,2,3) * rhostar * ((ci1*krho(imod_2)*phiga_tmp)*dapargads_tmp(i) &
                          & - (dphigads_tmp(i))*ci1*krho(imod_1)*aparga_tmp)) !zeta, s components
                    else
                      imod_2 = imod_1 - imod + 1
                      call dfielddx(APAR_GA_FIELD,imod_1,i,j,ispl,dapargadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      call dfieldds(APAR_GA_FIELD,imod_1,i,j,ispl,dapargads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      call dfielddx(PHI_GA_FIELD,imod_2,i,j,ispl,dphigadx_tmp) !dapargadx_xy_point_get(psi,zeta,i,j,ispl)
                      call dfieldds(PHI_GA_FIELD,imod_2,i,j,ispl,dphigads_tmp) !dapargads_xy_point_get(psi,zeta,i,j,ispl)
                      phiga_tmp = get_averaged_phi(imod_2,ix,i,j,ispl,fdisi)
                      aparga_tmp = get_averaged_apar(imod_1,ix,i,j,ispl,fdisi)
                      vchi_gradphi = vchi_gradphi &
                        & + 2.0 * vpgr(i,j,k,ispl) * vthrat(ispl) * ( &                        
                        & efun(ix,i,1,3) * rhostar * (conjg(dphigadx_tmp(ix))*dapargads_tmp(i) &
                          & - conjg(dphigads_tmp(i))*dapargadx_tmp(ix)) + & ! psi, s components
                        & efun(ix,i,2,3) * rhostar * (conjg(ci1*krho(imod_2)*phiga_tmp)*dapargads_tmp(i) &
                          & - conjg(dphigads_tmp(i))*ci1*krho(imod_1)*aparga_tmp)) !zeta, s components
                    end if
                    
                  end do
                  jdote_em_rhostar(k,j) = jdote_em_rhostar(k,j) &
                    & + signz(ispl) * real(conjg(vchi_gradphi) * dgrid(ispl) * fd) &
                    & * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)

                end if !nlapar.and.non_linear
              end if !not flux_tube
            end do; end do  ! Loop over s and imod (integration)
            jdote(k,j) = jdote_para_elec(k,j) + jdote_para_magnetic(k,j) + & 
              & jdote_perp_elec(k,j) + jdote_em(k,j) + &
              & jdote_perp_elec_rhostar(k,j) + jdote_em_rhostar(k,j)
          end do; end do ! Loop over vpar and mu
        else
          jdote(:,:) = 0.0
          jdote_para_elec(:,:) = 0.0
          jdote_para_magnetic(:,:) = 0.0
          jdote_perp_elec(:,:) = 0.0
          jdote_em(:,:) = 0.0
          jdote_perp_elec_rhostar(:,:) = 0.0
          jdote_em_rhostar(:,:) = 0.0
        end if

        !sum over s and x grid globally (values on other x processors are zero)
        call mpireduce_sum_inplace(jdote_para_elec, shape(jdote_para_elec), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_para_magnetic, shape(jdote_para_magnetic), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_perp_elec, shape(jdote_perp_elec), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_em, shape(jdote_em), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_perp_elec_rhostar, shape(jdote_perp_elec_rhostar), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote_em_rhostar, shape(jdote_em_rhostar), &
           & COMM_S_NE_X_NE)
        call mpireduce_sum_inplace(jdote, shape(jdote), &
           & COMM_S_NE_X_NE)

        luname="jdote_fsa_para_elec_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_para_elec,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 0)

        luname="jdote_fsa_para_magnetic_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_para_magnetic,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 1)

        luname="jdote_fsa_perp_elec_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_perp_elec,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 2)

        luname="jdote_fsa_em_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_em,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 3)

        luname="jdote_fsa_perp_elec_rhostar_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_perp_elec_rhostar,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 4)

        luname="jdote_fsa_em_rhostar_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote_em_rhostar,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 5)

        luname="jdote_fsa_pt"//trim(int2char_zeros(ipoint,2))// &
           & "_sp"//trim(int2char_zeros(ispg,2))// &
           & "_"//file_count_suffix
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & jdote,luname,ixg,1,ispg, binary_fmt, &
           & tag_range_start2 + (ipoint-1)*number_of_species*n_velocity_slices + (ispg-1)*n_velocity_slices + 6)
      end do
    end do

  end subroutine output_jdote_fs_pointwise


  !****************************************************************************
  !> return the real space value of the time-derivative
  !> of the gyro-averaged parallel magnetic
  !> potential at a specific psi zeta point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dapargadt_xy_point_get(psi,zeta,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use grid,             only : nmod, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use general,          only : gkw_abort
    use fields,           only : get_averaged_apar
    use dist,             only : fdisi
    use control,          only : flux_tube

    !> local parallel, mu, and species index
    integer, intent(in) :: i,j,is
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi,zeta
    !> the value of d<apar>/dt at the specified position is returned
    real :: dapargadt_xy_point_get

!    real                  :: ix_mid
    integer               :: ix, ixg, imod
    complex               :: dfield
    complex               :: aparga, daparga

    dfield= (0.0E0, 0.0E0)

    do imod = 1, nmod
      !ix_mid = real((n_x_grid+1)*0.5E0)
      !ix=nint(ix_mid+psi/lx*real(n_x_grid))
      if (flux_tube) then
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if

      ix = lrx(ixg)

      aparga = get_averaged_apar(imod,ix,i,j,is,fdisi)
      
      if (ixg < 1 .or. ixg > n_x_grid) then
        call gkw_abort('nonspectral: Error in dapargadt_xy_point_get located in diagnos_jdote.f90')
      end if

      daparga = (aparga - last_aparga(imod,ix,i,j,is)) / delta_time_jdote
      if (proc_subset(ixg,0,0,0,0)) then
        dfield = dfield + 2*daparga * &
           & exp(ci1*2*pi/lyn*(imod-1)*zeta)
      else ! psi point not on this processor
        dfield = 0.0
      end if
    end do

    dapargadt_xy_point_get = real(dfield)

  end function dapargadt_xy_point_get

  !****************************************************************************
  !> return the (complex) value of the time-derivative
  !> of the gyro-averaged parallel magnetic
  !> potential at a specific psi point.
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !>
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function dapargadt_x_point_get(psi,imod,i,j,is)
    use constants,        only : ci1, pi
    use mode,             only : lxn
    use grid,             only : n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use general,          only : gkw_abort
    use fields,           only : get_averaged_apar
    use dist,             only : fdisi
    use control,          only : flux_tube

    !> local parallel, mu, and species index
    integer, intent(in) :: i, j, is, imod
    !> radial and binormal position space coordinates (not indices!)
    real, intent(in) :: psi
    !> the value of d<apar>/dt at the specified position is returned
    complex :: dapargadt_x_point_get

    !real                  :: ix_mid
    integer               :: ix, ixg
    real                  :: dfield
    complex               :: aparga

    dfield=0.0E0

    !ix_mid = real((n_x_grid+1)*0.5E0)
    !ix=nint(ix_mid+psi/lx*real(n_x_grid))
    if (flux_tube) then
      ! We use the global number of radial points n_x_grid to determine the
      ! global ix for the chosen radial position
      ixg = nint(psi/lxn*real(n_x_grid))
    else
      ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
    end if
    ix = lrx(ixg)
    aparga = get_averaged_apar(imod,ix,i,j,is,fdisi)
    if (ixg < 1 .or. ixg > n_x_grid) then
      call gkw_abort('nonspectral: Error in dapargadt_xy_point_get')
    end if
    dapargadt_x_point_get = (aparga - last_aparga(imod,ix,i,j,is)) / delta_time_jdote
  end function dapargadt_x_point_get

end module diagnos_jdote
