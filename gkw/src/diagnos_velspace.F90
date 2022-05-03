!------------------------------------------------------------------------------
!> Output the distribution function or its moments as a 2D slice in 
!> velocity space.
!------------------------------------------------------------------------------
module diagnos_velspace

  implicit none 

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output
  public :: fdis_xy_point_get

  private

  logical, save, public :: lvelspace_output

  !>The number of points and then the x and zeta coordinates for the 
  !>velocity space outputs.
  integer, parameter :: nvpoints = 64
  integer, save, public :: npointsvel
  real, save, public :: psi_velspace(nvpoints)
  real, save, public :: zeta_velspace(nvpoints)

  !Outputs the full velocity space (in separate files) for 0,0 mode
  logical, save, public :: lfinvel

  integer, save :: i_momentum1 = -1, i_momentum2 = -1

  !> to avoid MPI dataraces due to gather operations inside
  !> gather_array inside velocity_slice_output, called in
  !> velspace_pointwise, tags are used. This number is the number of
  !> quantities output here in this routine:
  integer, parameter :: n_velocity_slices = 4
  !> the range of tags to be used by this diagnostic
  integer, save :: tag_range_start, tag_range_end_inkl
  integer, save :: tag_range_start2, tag_range_end_inkl2

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
      lvelspace_output = .false.
      
      npointsvel = 0
      psi_velspace(:) = 0.
      zeta_velspace(:) = 0.
      
      lfinvel = .false.

  end subroutine set_default_nml_values


  subroutine init(requirements)
    use mpiinterface, only : root_processor, register_tag_range
    use io, only : open_real_lu, ascii_fmt
    use grid, only : number_of_species
    logical, intent(inout) :: requirements(:,:)

    ! To keep the compiler quiet, as the array is not used so far (in this routine).
    if (.false.) write(*,*) requirements

    if(root_processor) then
      if(lvelspace_output) then


        call open_real_lu('Momentum1.dat', 'diagnostic/diagnos_velspace', (/ 3 /), &
           & ascii_fmt, i_momentum1)
        call open_real_lu('Momentum2.dat', 'diagnostic/diagnos_velspace', (/ 3 /), &
           & ascii_fmt, i_momentum2)
      end if
    end if

    if (npointsvel > 0) then
      call register_tag_range(npointsvel * number_of_species * n_velocity_slices, &
         & tag_range_start, tag_range_end_inkl)
      call register_tag_range(2, &
         & tag_range_start2, tag_range_end_inkl2)
    end if

  end subroutine init


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Routine for collision operator testing purposes
  !> fairly similar to velocity_space_output
  !----------------------------------------------------------------------------
  subroutine velspace(file_count)
    use io,             only : append_chunk, xy_fmt, ascii_fmt, binary_fmt
    use dist,           only : ifdis, fdisi
    use index_function, only : indx
    use grid,           only : nmu, nvpar, ns, nsp, number_of_species, gsp
    use control,        only : neoclassics
    use components,     only : adiabatic_electrons
    use geom,           only : bn, ints
    use mode,           only : ixzero, iyzero
    use velocitygrid,   only : intmu, intvp, vpgr, mugr
    use global,         only : int2char_zeros
    use matdat,         only : get_f_from_g
    use mpiinterface,   only : mpireduce_sum_inplace, root_processor
    use mpicomms,       only : COMM_SP_EQ_S_EQ

    use diagnos_generic, only : velocity_slice_output

    integer, intent(in) :: file_count
    integer :: i,j,k,ix 
    real :: dum, dumden, dum2, dumden2, dume, dume2
    real :: fdis
    character (len=14) :: luname
    real, allocatable, dimension(:,:) :: local_vpar_mu

    neocl: if (neoclassics) then

      ! select the first radial surface here 
      ix = 1 

      ! This works for a single mode only  
      !FJC_PARALLEL_X / NON-SPECTRAL - no obvious extension

      dum = 0.
      dumden = 0.
      dum2 = 0.
      dume = 0.
      dume2 = 0.

      do i=1,ns 
        do j=1,nmu
          do k=1,nvpar
            ! Zero mode should have no imaginary part, thus the conversion
            ! should be possible without loss of information.
            fdis = real(get_f_from_g(ixzero,iyzero,i,j,k,1,fdisi))
            dum=dum+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)*vpgr(i,j,k,1)*ints(i)*fdis
            dumden=dumden+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)*ints(i)*fdis
            dume=dume+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)*ints(i)* &
               & (vpgr(i,j,k,1)**2+2*mugr(j)*bn(ixzero,i))*fdis

            ! Zero mode should have no imaginary part, thus the conversion
            ! should be possible without loss of information.
            dum2=dum2+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)*vpgr(i,j,k,1)* & 
               & real(fdisi(indx(ifdis,1,1,i,j,k,2)))
            dumden2=dumden2+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)* & 
               & real(fdisi(indx(ifdis,1,1,i,j,k,2)))
            dume2=dume2+bn(ixzero,i)*intvp(i,j,k,1)*intmu(j)* & 
               & (vpgr(i,j,k,1)**2+2.E0*bn(ixzero,i)*mugr(j))* & 
               & real(fdisi(indx(ifdis,1,1,i,j,k,2)))
          end do
        end do
      end do

      call mpireduce_sum_inplace(dum,COMM_SP_EQ_S_EQ)
      call mpireduce_sum_inplace(dumden,COMM_SP_EQ_S_EQ)
      call mpireduce_sum_inplace(dume,COMM_SP_EQ_S_EQ)

      call mpireduce_sum_inplace(dum2,COMM_SP_EQ_S_EQ)
      call mpireduce_sum_inplace(dumden2,COMM_SP_EQ_S_EQ)
      call mpireduce_sum_inplace(dume2,COMM_SP_EQ_S_EQ)

      if(root_processor)then  
        write(*,*)'Momentum', dum, dumden, dume
        call append_chunk(i_momentum1, (/ dum,dumden,dume /), &
           & xy_fmt, ascii_fmt)
        call append_chunk(i_momentum2, (/ dum2,dumden2,dume2 /), &
           & xy_fmt, ascii_fmt)
      endif
    end if neocl


    ! allocate arrays to contain the full slice and local slice
    allocate(local_vpar_mu(nvpar,nmu))

    ! output the velocity space of species 1 at the point imod=1, ix=1, i=1
    i=1 
    do j=1,nmu
      do k=1,nvpar
        local_vpar_mu(k,j)=real(fdisi(indx(ifdis,1,1,i,j,k,1)))
      end do
    end do
    luname="IO"//trim(int2char_zeros(file_count,8))//".dat"
    call velocity_slice_output('diagnostic/diagnos_velspace', &
       & local_vpar_mu,luname,1,1,1, binary_fmt, tag_range_start2)

    ! note: The following code assumes that the electrons are last, i.e. have
    !       highest index.
    if(.not. adiabatic_electrons) then
      do j=1,nmu
        do k=1,nvpar
          local_vpar_mu(k,j)=real(fdisi(indx(ifdis,1,1,i,j,k,nsp)))
        end do
      end do
      luname="EL"//trim(int2char_zeros(file_count,8))//".dat"
      call velocity_slice_output('diagnostic/diagnos_velspace', &
         & local_vpar_mu,luname,1,1,number_of_species, binary_fmt, tag_range_start2+1)
    end if

    ! deallocate the temporary arrays
    if (allocated(local_vpar_mu)) deallocate(local_vpar_mu)

  end subroutine velspace


  !****************************************************************************
  !> Routine to output velocity space at selected (x,y) points.
  !> These points are selected using psi_velspace and zeta_velspace
  !> and npointsvel in the diagnostic namelist.  The values can go from 
  !> psi = -/+ LX/2 and zeta = -/+ LY/2 relative to the xphi and yphi grids
  !> where 0,0 is centred on the island center.
  !> The values are integrated along the parallel direction (can be changed)
  !> and the integrand (currently radial / poloidal drifts) can also be changed.
  !> Call from subroutine diagnostic_naverage
  !----------------------------------------------------------------------------
  subroutine velspace_pointwise(file_count)
    use mpiinterface,    only : mpiallreduce_sum, mpireduce_sum_inplace
    use mpiinterface,    only : mpibarrier
    use mpicomms,        only : COMM_S_NE
    use control,         only : time
    use geom,            only : ints, bn, dfun, isg_lfs
    use velocitygrid,    only : intvp, intmu, mugr, vpgr
    use global,          only : int2char_zeros
    use grid,            only : ns, nmu, nvpar, n_x_grid, gs
    use grid,            only : lsp, number_of_species, lrx, proc_subset
    use components,      only : signz, tmp, tearingmode
!    use dist,            only : f_EP
    use dist,            only : fmaxwl
    use diagnos_generic, only : velocity_slice_output
    use io,              only : binary_fmt
    use mode,            only : lxn, lyn
    use tearingmodes,    only : omega_rot
    use constants,       only : pi

    integer, intent(in) :: file_count
    real :: psi, zeta, fd, vdriftr, vdriftp
    integer :: i,j,k,ispl,ispg, ipoint, ixg, ix
    character (len=64) :: luname_rad
    character (len=64) :: luname_pol
    character (len=64) :: luname_eq
    character (len=64) :: luname_vel_lfs

    !>two different vgrid_slices because, for example, if
    !>you want to use the pointwise diagnostic for the fluxes, you have two
    !>components to consider
    real, allocatable, save, dimension(:,:) :: vgrid_slice_rad, vgrid_slice_pol, vel_lfs
    !>Equilibrium distribution
    real, allocatable, save, dimension(:,:) :: Feq
    
    character (len=6), save :: file_count_suffix
    file_count_suffix=trim(int2char_zeros(file_count,6))

    ! todo for global:
    !   * determine radial coordinate index correctly from psi. Note that box
    !     goes from psil to psih
    !   * use also the dfieldds subroutine to calculate the energy exchanges
    !     between particles and waves in the parallel direction. This
    !     contribution is O(rho*) for flux tube runs
    !   * ...likely other things
    
    do ipoint=1,npointsvel ! Loop over the output points you are interested in

      !allocate array to contain the local slice
      if (.not. allocated(vgrid_slice_rad))    allocate(vgrid_slice_rad(nvpar,nmu))
      if (.not. allocated(vgrid_slice_pol))   allocate(vgrid_slice_pol(nvpar,nmu))
      if (.not. allocated(Feq)) allocate(Feq(nvpar,nmu))
      if (.not. allocated(vel_lfs)) allocate(vel_lfs(nvpar,nmu))

      !The coordinates of the points you are interested in have to
      !be introduced referring to the xphi-yphi grid AFTER the
      !reorganization of the matrix (see below)
      psi=psi_velspace(ipoint)
      zeta=zeta_velspace(ipoint)

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
      ! Also, psi is the radial position, normalised to rho_ref.
      ! This is to place psi=0 at the centre of the box:
      psi=psi+lxn/2
      ! We use the global number of radial points n_x_grid to determine the
      ! global radial index for the chosen radial position
      ixg = nint(psi/lxn*real(n_x_grid))
      ix = lrx(ixg)

      do ispg = 1, number_of_species ! Beginning of loop over global species
        ! We need to determine the local species number
        ! (in the sense of processor)
        ispl=lsp(ispg)

        !if that species is not on this processor skip to next species
        if (proc_subset(ixg,0,0,0,ispg)) then
          vgrid_slice_rad = 0.0
          vgrid_slice_pol = 0.0
          Feq = 0.0
          vel_lfs = 0.0
          
          do i=1,ns ! Beginning of loop over local parallel coordinate
            ! After this loop we need to communicate between processors
            ! to complete the flux surface average
            do j = 1, nmu; do k=1,nvpar ! Loop over velocity points

              ! But can also be re-calculated as follows,
              ! without the beta contribution
              vdriftr = (vpgr(i,j,k,ispl)**2+(mugr(j)*bn(ix,i)))* &
                 & dfun(ix,i,1)*tmp(ix,ispl)/signz(ispl)
              vdriftp = (vpgr(i,j,k,ispl)**2+(mugr(j)*bn(ix,i)))* &
                 & dfun(ix,i,2)*tmp(ix,ispl)/signz(ispl)


              ! This diagnostic refers to the quantity vd*f in the
              ! velocity space. If you want for example the density,
              ! put vdriftr=vdriftp=1.

              ! Get the real value of fdis at the point (zero if not on this x-processor)
              !if (energetic_particles) then
              
              fd = fdis_xy_point_get(psi,zeta,i,j,k,ispl)
              
              !else
              !  fd = 1.E0
              !end if

              !do whatever else you need to it
              !Element of density = distribution*velocity space Jacobian
              vgrid_slice_rad(k,j) = vgrid_slice_rad(k,j) + vdriftr*fd* &
                 & bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
              vgrid_slice_pol(k,j) = vgrid_slice_pol(k,j) + vdriftp*fd* &
                 & bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
              if(proc_subset(ixg,0,0,0,0)) then
                Feq(k,j) = Feq(k,j) + fmaxwl(ix,i,j,k,ispl)* &
                   & bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
                !if (energetic_particles) then
                !  Feq(k,j) = Feq(k,j) + f_EP(ix,i,j,k,ispl)*ints(i)
                !end if
              end if
              
              !velocity space slice of fdisi at lfs (s=0)
              if(gs(i) == isg_lfs) then
                vel_lfs(k,j) = fd * bn(ix,i)*intvp(i,j,k,ispl)*intmu(j)*ints(i)
              end if
              
            end do; end do
          end do

          call mpireduce_sum_inplace(vgrid_slice_rad, shape(vgrid_slice_rad), &
             & COMM_S_NE)
          call mpireduce_sum_inplace(vgrid_slice_pol, shape(vgrid_slice_pol), &
             & COMM_S_NE)
          call mpireduce_sum_inplace(Feq,shape(Feq),COMM_S_NE)

        end if


        !construct luname based on your point number and species number 
        luname_rad="vd_rad.pt"//trim(int2char_zeros(ipoint,2))// &
           & ".sp"//trim(int2char_zeros(ispg,2))// &
           & ".dat"//file_count_suffix

        luname_pol="vd_pol.pt"//trim(int2char_zeros(ipoint,2))// &
           & ".sp"//trim(int2char_zeros(ispg,2))// &
           & ".dat"//file_count_suffix

        luname_eq="Feq.pt"//trim(int2char_zeros(ipoint,2))// &
           & ".sp"//trim(int2char_zeros(ispg,2))//&
           & ".dat"//file_count_suffix

        luname_vel_lfs="vel_lfs.pt"//trim(int2char_zeros(ipoint,2))// &
           & ".sp"//trim(int2char_zeros(ispg,2))// &
           & ".dat"//file_count_suffix

        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & vgrid_slice_rad,luname_rad,ixg,1,ispg, binary_fmt, &
           & tag_range_start - 1 + (ipoint-1)*number_of_species*n_velocity_slices + &
           & (ispg-1)*n_velocity_slices + 1)
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & vgrid_slice_pol,luname_pol,ixg,1,ispg, binary_fmt, &
           & tag_range_start - 1 + (ipoint-1)*number_of_species*n_velocity_slices + &
           & (ispg-1)*n_velocity_slices + 2)
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & Feq,luname_eq,ixg,1,ispg, binary_fmt, &
           & tag_range_start - 1 + (ipoint-1)*number_of_species*n_velocity_slices + &
           & (ispg-1)*n_velocity_slices + 3)
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & vel_lfs,luname_vel_lfs,ixg,isg_lfs,ispg, binary_fmt, &
           & tag_range_start - 1 + (ipoint-1)*number_of_species*n_velocity_slices + &
           & (ispg-1)*n_velocity_slices + 4)
      end do !species
    end do ! your points

  end subroutine velspace_pointwise

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Output velocity grids at s point ipar_global
  !> Output velocity space diagnostic for all species at s point ipar_global
  !> Assumes grids already written by velocity_space_output
  !> At present only useful for mode_box false runs - plots mode(1,1)
  !----------------------------------------------------------------------------
  subroutine final_velocity_output

    use global,         only : int2char_zeros
    use dist,           only : ifdis, fdisi
    use index_function, only : indx
    use grid,           only : nmu,nvpar, lsp, ls, number_of_species 
    use grid,           only : proc_subset, n_s_grid
    use mode,           only : ixzero
    use diagnos_generic, only : velocity_slice_output
    use io, only : binary_fmt
    use mpiinterface, only : register_tag_range

    integer :: isg,j,k, ispg, ispl
    real, allocatable, dimension(:,:) :: local_vpar_mu
    character (len=18) ::  luname
    integer :: tag_range_start3, tag_range_end_inkl3
    
    call register_tag_range(number_of_species*n_s_grid, &
       & tag_range_start3, tag_range_end_inkl3)

    ! allocate arrays to contain the full slice and local slice
    allocate(local_vpar_mu(nvpar,nmu))

    !Loop over all species
    species: do ispg=1,number_of_species
      !local species
      ispl=lsp(ispg)

      !Loop over s points
      spoints: do isg=1,n_s_grid

        if(proc_subset(0,isg,0,0,ispg)) then
          ! local processor real part of fdisi
          do j=1,nmu
            do k=1,nvpar
#define FVELSPACE_IMOD 1
              local_vpar_mu(k,j)=real(fdisi(indx(ifdis,FVELSPACE_IMOD, &
                 & ixzero,ls(isg),j,k,ispl)))
            end do
          end do
        end if

        luname="vs.s"//trim(int2char_zeros(isg,3))// &
           & ".sp"//trim(int2char_zeros(ispg,2))//".dat"
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & local_vpar_mu,luname,ixzero,isg,ispg, binary_fmt,&
           & tag_range_start3+n_s_grid*(ispg-1) + isg-1)

      end do spoints
    end do species

    ! deallocate the temporary arrays
    if (allocated(local_vpar_mu)) deallocate(local_vpar_mu)

  end subroutine final_velocity_output

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> Output velocity grids at x point ipar_global
  !> Output velocity space diagnostic for all species at s point ipar_global
  !> Assumes grids already written by velocity_space_output
  !> At present only useful for mode_box false runs - plots mode(1,1)
  !----------------------------------------------------------------------------
  subroutine final_velocity_output2

    use global,         only : int2char_zeros
    use dist,           only : ifdis, fdisi
    use index_function, only : indx
    use grid,           only : n_x_grid, nmu,nvpar, lsp, lrx, number_of_species
    use grid,           only : proc_subset, n_x_grid, ls
    use diagnos_generic, only : velocity_slice_output
    use io, only : binary_fmt
    use mpiinterface, only : register_tag_range

    integer :: ixg,j,k, ispg, ispl, isg
    real, allocatable, dimension(:,:) :: local_vpar_mu
    character (len=18) ::  filename
    integer :: tag_range_start4, tag_range_end_inkl4

    call register_tag_range(number_of_species*n_x_grid, &
       & tag_range_start4, tag_range_end_inkl4)
    
    ! allocate arrays to contain the full slice and local slice
    allocate(local_vpar_mu(nvpar,nmu))

    isg = 1
    
    !Loop over all species
    species: do ispg=1,number_of_species
      !local species
      ispl = lsp(ispg)

      !Loop over x points
      xpoints: do ixg=1,n_x_grid

        if(proc_subset(ixg,isg,0,0,ispg)) then
          ! local processor real part of fdisi
          do j=1,nmu
            do k=1,nvpar
#define FVELSPACE_IMOD 1
              local_vpar_mu(k,j)=real(fdisi(indx(ifdis,FVELSPACE_IMOD, &
                 & lrx(ixg),ls(isg),j,k,ispl)))
            end do
          end do
        end if

        filename="vs.x"//trim(int2char_zeros(ixg,3))// &
           & ".sp"//trim(int2char_zeros(ispg,2))//".dat"
        call velocity_slice_output('diagnostic/diagnos_velspace', &
           & local_vpar_mu,filename,ixg,isg,ispg, binary_fmt, &
           & tag_range_start4 + n_x_grid*(ispg-1) + ixg-1)

      end do xpoints
    end do species

    ! deallocate the temporary arrays
    if (allocated(local_vpar_mu)) deallocate(local_vpar_mu)

  end subroutine final_velocity_output2


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(lvelspace_output,1)
    call mpibcast(npointsvel,1)
    call mpibcast(psi_velspace,nvpoints)
    call mpibcast(zeta_velspace,nvpoints)
    call mpibcast(lfinvel,1)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use control, only : flux_tube
    use mode, only : mode_box, lxn, lyn
    use general, only : gkw_warn, gkw_abort
    integer :: ipoint
    real :: psi, zeta
    
    if(.not.lvelspace_output) return

    if (.not.mode_box) then
      call gkw_warn('the lvelspace_output diagnostic needs to do a Fourier trafo&
         &, but this is not possible for mode_box=F. Switched off')
      ! the trafo is in fdis_xy_point_get
      lvelspace_output = .false.
    end if

    do ipoint=1,npointsvel
      psi=psi_velspace(ipoint)
      if (psi > lxn/2. .or. psi < -lxn/2.) then
        write(*,*) 'psi_velspace:', psi
        call gkw_abort('bad psi_velspace point selection')
      end if

      zeta=zeta_velspace(ipoint)
      if (zeta > lyn/2. .or. zeta < -lyn/2.) then
        write(*,*) 'zeta_velspace:', zeta, 'not in range', -lyn/2., lyn/2.
        call gkw_abort('bad zeta_velspace point selection')
      end if
    end do

    if (.not.flux_tube) then 
      call gkw_warn('Called velspace_pointwise with flux_tube = .false.,&
         & The correct dependence of T on x has never been implemented')
      lvelspace_output = .false.
    endif

  end subroutine check


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    !if (allocated(vgrid_slice)) deallocate(vgrid_slice)
    !if (allocated(vgrid_slice2)) deallocate(vgrid_slice2)
    !if (allocated(vgrid_slice_av)) deallocate(vgrid_slice_av)
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
    use components, only : tm_drive

    ! velocity-space diagnostic #3:
    if(lfinvel) call final_velocity_output
    if(lfinvel.and.tm_drive) call final_velocity_output2

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output(file_count)
    use mode, only : mode_box
    use diagnos_generic, only : xy_estep
    integer, intent(in) :: file_count

    ! diagnostic #1
    !Collision operator testing routines.
    if (lvelspace_output) then
      write(*,*)'Calculate momentum'
      call velspace(file_count)
    end if

    ! if it is not modebox=T then the fourier trafo in
    ! fdis_xy_point_get does not make sense.
    if (.not.mode_box) return 

    ! diagnostic #2
    if (xy_estep) then
      if ((npointsvel > 0).and.(lvelspace_output)) call velspace_pointwise(file_count)
    end if
  end subroutine output

    
  !****************************************************************************
  !> return the real space value of f at a specific psi zeta point 
  !> the value of f array for vpar, mu, s, sp, zeta, psi is returned
  !> SLOW! not to be used for loops over many x points - use FFTs instead!
  !> The points i,j,k,is are local
  !> If the psi point is not on the local processor (e.g. parallel_x)
  !> then zero is returned (i.e. an allreduce / selection must be performed
  !> outside this routine)
  !----------------------------------------------------------------------------
  function fdis_xy_point_get(psi,zeta,i,j,k,is)
    use dist,             only : fdisi
    use constants,        only : ci1, pi
    use mode,             only : lxn, lyn
    use control,          only : spectral_radius, flux_tube
    use grid,             only : nx, nmod, lx, n_x_grid, lrx, proc_subset
    use grid,             only : psil, psih
    use matdat,           only : get_f_from_g
    use general,          only : gkw_abort

    integer, intent(in) :: i,j,k,is
    real, intent(in) :: psi,zeta

    real    :: fdis_xy_point_get
!    reak    :: ix_mid
    integer :: ix, ixg, imod, helpint
    complex :: fdis, dum 

    fdis = (0.0,0.0)

    if (spectral_radius) then    !Manually do the fourier sum in 2D
      do imod = 1, nmod
        do ix = 1, nx  
          ! dum is the distribution without A_par contribution  
          dum = get_f_from_g(imod,ix,i,j,k,is,fdisi)    

          ! This remaps the island to 0,0
          helpint=ix-1-(nx-1)/2

          ! should have used krho and kxrh here
          fdis = fdis + 2*dum*exp(ci1*2*pi/lyn*(imod-1)*zeta + &
             & ci1*2*pi*helpint/lx*psi)
        end do
      end do

    else
      ! Manually do the fourier sum in 1D, choose the nearest point in x .
      ! Note that island is at centre in nonspectral.
      !ix_mid = real((n_x_grid+1)*0.5E0)      
      !ixg=nint(ix_mid+psi/lx*real(n_x_grid))
      if (flux_tube) then
        ! We use the global number of radial points n_x_grid to determine the
        ! global ix for the chosen radial position
        ixg = nint(psi/lxn*real(n_x_grid))
      else
        ixg = nint((psi-psil)*real(n_x_grid)/(psih-psil)+0.5E0)
      end if

      if (ixg < 1 .or. ixg > n_x_grid) then
        call gkw_abort('nonspectral: Error in fdis_xy_point_get located in diagnos_velspace.F90')
      end if

      if (proc_subset(ixg,0,0,0,0)) then
        do imod = 1, nmod 
          ! dum is the distribution without A_par contribution  
          dum = get_f_from_g(imod,lrx(ixg),i,j,k,is,fdisi)    

          ! should have used krho here
          fdis = fdis + 2*dum*exp(ci1*2*pi/lyn*(imod-1)*zeta)
        end do
      else ! psi point not on this processor    
        fdis = (0.0,0.0)
      end if
    end if

    fdis_xy_point_get = real(fdis) 

  end function fdis_xy_point_get


end module 
