!-----------------------------------------------------------------------------
!> Contains initialisation of the magnetic island structure.
!> The island structure is initialised here but applied in fields
!-----------------------------------------------------------------------------
module tearingmodes

  implicit none

  private
  
  public :: tearingmodes_init
  
  public :: initialise_island,islandstruct, isl_phi_indx
  public :: islandindx, omega_rot, isl_phi, imodisland
  
  !The poloidal mode where the island is
  integer, save :: imodisland

  !Normalised island rotation frequency
  real, save :: omega_rot

  !Island fields structures and indexing arrays
  complex, save, allocatable :: islandstruct(:,:)
  complex, save, allocatable :: isl_phi(:,:,:)
  integer, save, allocatable :: islandindx(:,:)
  integer, save, allocatable :: isl_phi_indx(:,:,:)
  
contains

  !---------------------------------------------------------------------------
  !> 
  !---------------------------------------------------------------------------
  subroutine tearingmodes_init

    ! The magnetic island is initialised as a perturbation in the
    ! parallel vector potential.
    call initialise_island

  end subroutine tearingmodes_init
  
!-----------------------------------------------------------------------------
!> Subroutine that initialises the parallel vector potential 
!> to have a magnetic island structure.
!> The island structure is initialised here but applied in fields.
!-----------------------------------------------------------------------------

subroutine initialise_island
  use general,          only : gkw_abort, gkw_warn
  use mpiinterface,     only : root_processor
  use grid,             only : nmod,nx,ns, lx, n_x_grid, gx
  use dist,             only : iapar
  use index_function,   only : indx
  use geom,             only : sgr, dxgr, q, bmin       
  use geom,             only : eps, shift_end_grid, shat
  use mode,             only : mode_box, ixzero, ikxspace, krho
  use mode,             only : lshat_zero
  use components,       only : wstar, tearingmode, isl_shear
  use components,       only : isl_Ls, isl_rot_freq, isl_mode
  use components,       only : psi_0, delta_psi_0    
  use constants,        only : pi
  use control,          only : spectral_radius, flux_tube
    
  ! integers for the loop over all grid points 
  integer :: i, p, ierr, ikxspace_local = 0
  integer :: ix

  ! Dummy variables 
  real    :: balloonpos
  complex :: dumbuf
  real :: psi_tmp, ix_mid

  ! Normalised amplitude of apar for the magnetic island
  real :: aparamp 

  !The maximum number of harmonics in the sum to create the
  !magnetic island
  integer :: iplusmode,iminusmode, en = 0

  !The poloidal mode where the magnetic island perturbation is
  !placed. Either the first or second.  First if only one mode
  !the second if greater (when mode_box=true) because the first
  !mode is alway zero frequency.
  integer :: modnum
  
  ! Two parameters are introduced to force the A|| perturbation to zero at the
  ! edge of simulation domain (in the radial direction). Be careful that the 
  ! radial discretization is good enough.
  ! THESE SHOULD BE DOCUMENTED IN THE INPUT FILE AND IN COMPONENTS 

  omega_rot = 0
  
  if (.not.tearingmode) return
  if (spectral_radius) then
    if (mode_box) then
      ikxspace_local  = ikxspace          
    else !if not mode_box ikxspace has no meaning in the rest of the code
      ikxspace_local  = 1
    end if
    !Comment since it is useless
    !islandtotpoints = nx*ns
  end if       
    
  ! Allocate for the island structure and the island index.
  allocate(islandstruct(nx,ns),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate the array &
                         & islandstruct in island initialisation')
  allocate(islandindx(nx,ns),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate the array &
                             & islandindx in island initialisation')

  ! Magnetic island perturbation added to apar
  ! The whole following section just initialises the parallel vector potential
  ! The distribution function itself isnt initialised until later.

  ! Mode amplitudes for the magnetic island
  ! There should be a further factor of Bn/Rn here but both these are set
  ! to one (in s-alpha geom only) and are therefore neglected.
  ! aparamp = 0.25E0*bmin*wstar**2/(q*pi)y

  if(.not. flux_tube) then
    call gkw_warn('bmin radial dependency has been ignored concerning tearing &
       &modes: bmin(1)')
  end if
  if (lshat_zero) then
    aparamp = 0.5E0*0.25E0*bmin(1)*isl_shear*wstar**2
  else
    if (spectral_radius) then
      aparamp = ikxspace_local*0.25E0*bmin(1)*eps*wstar**2/(q*q*pi)
    else
      aparamp = 0.25E0*bmin(1)*shat*wstar**2/q
    end if
  end if

  omega_rot = isl_rot_freq

  ! Writing useful information
  if(root_processor) then
    write(*,*)'Tearing mode', tearingmode
    if (isl_rot_freq.gt.1e-10) then
      write(*,*)'Island rotation',  isl_rot_freq
    end if
    write(*,*)'Magnetic island structure added'
    write(*,*)'Mode damping length set to', isl_Ls
    write(*,*)'Amplitude of island perturbation is',aparamp    
    ! The two methods below (nmod.eq.1) and else are 
    ! equivalent due to the ballooning transform
  end if
    
  if (nmod.eq.1) then
    imodisland=1
    en = 1
  else if (nmod.gt.1) then       
    imodisland=isl_mode
    en = imodisland-1
    ! FJC: why is this commented ?
    ! DZ : I think the amplitude should not be divided by n
    ! aparamp = aparamp/en
    ! THE MANUAL SHOULD BE CORRECTED
    if (root_processor) then
      write(*,*)'There are ', imodisland-1, ' islands in the box' 
    end if       
  end if    
  if (imodisland.gt.nmod) call gkw_abort('Tearing mode: The poloidal mode &
     & of the island is defined as must be less than nmod')

  if (spectral_radius) then
       
    ! The constant psi approximation is broken so as to reduce the effect of
    ! discontinuities on the boundary.  Here the width of the gaussian is 
    ! defined to be approximately 8 radial modes.    
    do i=1,ns
      if (lshat_zero) then 
        islandstruct(ixzero,i) = aparamp*(1.E0,0.E0)                  
      else
        if (abs(sgr(i)).lt.1e-10) then
          dumbuf = pi*aparamp*(1.E0,0.E0)
          ! Set island amplitude in the Apar field for ixzero mode
          ! And in the island struct that will be reapplied
          islandstruct(ixzero,i)= dumbuf
        else             
          balloonpos = ikxspace_local*en*sgr(i)
          dumbuf = aparamp*exp(-(balloonpos/isl_Ls)**2)*sin(pi*balloonpos) &
                 &   *(1.E0,0.E0)/balloonpos
          islandstruct(ixzero,i)=dumbuf
          islandindx(ixzero,i)= indx(iapar,imodisland,ixzero,i)
        end if
      end if
      islandindx(ixzero,i)= indx(iapar,imodisland,ixzero,i)
    end do
         
    ! number of non-zero positive frequency radial modes
    modnum = (nx-1)/2
       
    if (nx.gt.1) then
      ! This is loop over x modes
      ! Initialise symmetric in kx.
      if (lshat_zero) then
        do p = 1,modnum
          iplusmode = ixzero+p
          iminusmode = ixzero-p
          do i=1,ns
            islandstruct(iplusmode,i) = 0.             
            islandstruct(iminusmode,i) = 0.             
            islandindx(iminusmode,i)= indx(iapar,imodisland,iminusmode,i)
            islandindx(iplusmode,i)=  indx(iapar,imodisland,iplusmode,i)
          end do
        end do
      else
        do p = 1,modnum
          iplusmode = ixzero+p
          iminusmode = ixzero-p
          do i=1,ns
            balloonpos = (ikxspace_local*en*sgr(i)+1.E0*p)
            dumbuf = (1.E0,0.E0)*aparamp*exp(-(balloonpos/isl_Ls)**2)* &
                   & sin(pi*balloonpos) /balloonpos                   
            ! Initialise island amplitude in connected x modes               
            !in the island struct that will be re-applied
            islandstruct(iplusmode,i) = dumbuf
            islandindx(iplusmode,i)= indx(iapar,imodisland,iplusmode,i)
                   
            balloonpos = (ikxspace_local*en*sgr(i)-1.E0*p)
            dumbuf = (1.E0,0.E0)*aparamp*exp(-(balloonpos/isl_Ls)**2)* &
                   & sin(pi*balloonpos) /balloonpos    
            ! And do the same in the kx refelction
            ! and in the island struct that will be re-applied
            islandstruct(iminusmode,i) = dumbuf
            islandindx(iminusmode,i)= indx(iapar,imodisland,iminusmode,i)
          end do
        end do
      end if
    end if
       
  else  ! nonspectral case

    ix_mid = real((n_x_grid+1)*0.5E0)
    if (lshat_zero) then          
      do i=1,ns; do ix=1,nx             
        psi_tmp = lx*(real(gx(ix))-ix_mid)/n_x_grid
        islandstruct(ix,i) = aparamp*(1.E0,0.E0)                      &
               &     *0.25E0*(1+tanh((psi_tmp+psi_0)/delta_psi_0))    &
               &     *(1-tanh((psi_tmp-psi_0)/delta_psi_0))
        islandindx(ix,i)=  indx(iapar,imodisland,ix,i)                
      end do; end do
    else
      do i=1,ns; do ix = 1,nx
        psi_tmp = dxgr * (real(gx(ix))-ix_mid)
        dumbuf = aparamp*((1.E0,0.E0)*cos(krho(imodisland)*           &
               &     shift_end_grid(ix)*sgr(i)) - (0.E0,1.E0)*        &
               &     sin(krho(imodisland)*shift_end_grid(ix)*sgr(i))) &
               &     *0.25E0*(1+tanh((psi_tmp+psi_0)/delta_psi_0))    &
               &     *(1-tanh((psi_tmp-psi_0)/delta_psi_0))
        islandstruct(ix,i)= dumbuf
        islandindx(ix,i)  = indx(iapar,imodisland,ix,i)
      end do; end do
    end if
  end if
    
  call init_tear_zero_epar
    
end subroutine initialise_island
    
  !---------------------------------------------------------------------------
  !> Initalise the zero_epar correction for rotating islands:
  !> We want to impose the analytic prediction for the
  !> electrostatic potential (i.e. impose Epar=0)
  !> This is to match analytics, but is not physical.
  !> Spectral only.
  !---------------------------------------------------------------------------
  subroutine init_tear_zero_epar
    use constants, only : pi, ci1
    use general, only : gkw_abort
    use control, only : spectral_radius
    use grid, only : ns, nmod, n_x_grid, lx, n_y_grid
    use mode, only : krho, kxrh
    use geom, only : sgr
    use dist, only : iphi
    use mode, only : lyinv, lyn !Note lyn /= 1/lyinv 
    use components, only : tear_zero_epar
    use fft, only : four2d_real, working_fft_library, FFT_FORWARD
    use grid, only : jinv_flexible
    use components, only : wstar
    use index_function,   only : indx
!    use io, only : xy_fmt, ascii_fmt
!    use mpiinterface, only : root_processor

    integer :: i, k, j, ierr
    real :: xr, yr
    real :: omega_pertlab, h_omega
    complex, allocatable :: phic(:,:)
    real, allocatable :: phir(:,:)
    integer :: jinv_flex
    integer :: n_x_grid_highres
    integer :: n_y_grid_highres, nmod_highres
       
    if (.not. tear_zero_epar) return
    if (.not.spectral_radius) then
      call gkw_abort('tear_zero_epar correction is currently only implemented &
         & with spectral')
    end if
    if (.not. working_fft_library) then
      call gkw_abort('tear_zero_epar requires an fft library')
    end if

    ! as one imposes a function with discontinuity in position space,
    ! a high resolution may help to reduce ringing artifacts.  One
    ! could choose an arbitrary, sufficiently large factor between the
    ! normal gridsizes and *_highres, for example 2.0.
    ! In order not to break the impose_tm2_nl testcase, choose
    n_x_grid_highres = int(1.5610 * n_x_grid)
    n_y_grid_highres = int(1.8 * n_y_grid)
    
    nmod_highres = floor(n_y_grid_highres/2.0) + 1
    
    allocate(phic(nmod_highres,n_x_grid_highres),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate the array &
       & phic in island initialisation')
    allocate(phir(n_y_grid_highres,n_x_grid_highres),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate the array &
       & phir in island initialisation')

    allocate(isl_phi(nmod,n_x_grid,ns),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate the array &
       & isl_phi in island initialisation')
    allocate(isl_phi_indx(nmod,n_x_grid,ns),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate the array &
       & isl_phi_indx in island initialisation') 

    ! fill perpendicular position space slice of the potential, called phir
    do i=1,ns

      do k=1,n_y_grid_highres; do j=1,n_x_grid_highres

        ! box coordinate on scale [-0.5,0.5]
        ! Island centre is at -0.5,-0.5 in these coordinates?
        ! Also phase shift in kx/ky after the FFT to make a translation
        xr=real(j-(n_x_grid_highres)/2)/ n_x_grid_highres
        yr=real(k-n_y_grid_highres/2)/ n_y_grid_highres

        ! the function omega_pertlab represents the dimensionless
        ! flux-surface label in presence of a magnetic island, 
        ! typically indicated with capital Omega
        omega_pertlab=2*xr*xr*lx**2/wstar**2-cos(2*pi*(yr+0.5-xr*sgr(i)))

        if (omega_pertlab > 1 ) then
          h_omega = SIGN(1.0,xr)*wstar/SQRT(2.0)/lx*(SQRT(omega_pertlab)-1)
          ! h_omega is the well-known profile function h(\Omega)
          ! Here we employ the simple form suggested by Smolyakov

          phir(k,j)= omega_rot*lx*lyn/pi*(xr-h_omega)*EXP(-(omega_pertlab-1)/20)
          ! Artificial exp smoothing to reduce effects at the edge =>
          ! The condition E// = 0 might not be satisfied far away from the island
        else
          phir(k,j)= omega_rot*lx*lyn/pi*xr
        end if

        ! Potential problems when ixzero not defined?

      end do; end do

      ! Write a file to check the form of your function at s=0
      ! if (abs(sgr(i))<1e-10 .and. root_processor) then
      !   call output_array('phir', 'diagnostic/diagnos_tearing', phir, 'F', &
      !     & xy_fmt, ascii_fmt)
      ! end if

      call four2d_real(phir,phic,FFT_FORWARD)

      ! Normalise
      phic=phic/(n_x_grid_highres*n_y_grid_highres)

      ! if (abs(sgr(i))<1e-10 .and. root_processor) then
      !   call output_array('jinv_flex', &
      !     & 'diagnostic/diagnos_tearing', &
      !     & (/ (1.0*jinv_flexible(n_x_grid_highres, n_x_grid_highres,j), &
      !     & j=1,n_x_grid_highres) /), 'F', &
      !     & xy_fmt, ascii_fmt)
      ! end if

      ! Put into island phi array
      do k=1,nmod
        do j=1,n_x_grid_highres
          jinv_flex = jinv_flexible(n_x_grid,n_x_grid_highres,j)
          ! if (jinv(j) /= ixzero .and. jinv(j).ne.0) then 
          if (jinv_flex.ne.0) then
            ! Set phi to the value prescribed, and translate in LX/LY
            ! box (depends on coords used above)
            isl_phi(k,jinv_flex,i)=phic(k,j) * &
               & exp(0.5*ci1*(lx*kxrh(jinv_flex)+krho(k)/lyinv))
            isl_phi_indx(k,jinv_flex,i)=indx(iphi,k,jinv_flex,i)
          end if
        end do
      end do

      ! Transform back, simply to check if you wish
      ! This will break the array in phic, uses it as workspace
      ! Write a file to check the form of your function at s=0
      ! if (abs(sgr(i))<1e-10 .and. root_processor) then
      !   call four2d_real(phir,phic,FFT_INVERSE)
      !   call output_array('phir_backtransformed', &
      !     & 'diagnostic/diagnos_tearing', phir, 'F', &
      !     & xy_fmt, ascii_fmt)
      ! end if
      
    end do !s

    deallocate(phir)
    deallocate(phic)
    
  end subroutine init_tear_zero_epar 
        

end module tearingmodes
