!-----------------------------------------------------------------------------
!> Calculates the nonlinear term III with FFT pseudospectral method.
!> For nonlinear runs, the code spends 40% (or much more) of time here.
!> This module is optimised, contains OpenMP, and (optional) perflib timings
!-----------------------------------------------------------------------------
module non_linear_terms

  implicit none

  private

  public :: add_non_linear_terms, jind, jinv, mphi, mphiw3, mrad
  public :: nonlinear_allocate, nonlinear_init, nl_initialised
  public :: entropy_radial
  public :: get_extended_firstdim_fft_size
  public :: mbd, ibd
  
  !> Index array for the storage of the kx modes in the arrays for the FFT
  !> jind(nx)
  integer, save, allocatable :: jind(:)

  !> Inverse index array for the storage of the kx modes in the arrays for the
  !> FFT. jinv(mrad)
  integer, save, allocatable :: jinv(:)

  !>Nested loop indexing array
  integer, save, allocatable :: i3loop(:,:)
  ! arrays for the FFT

  ! Throughout the code:
  ! * x is used to refer to the radial (psi) direction
  ! * y is used to refer to the perpendicular (zeta) direction within the flux
  !     surface, This direction can be called 'polodial' but is also referred
  !     to as 'toroidal'
  ! * Gradient d/dx in k space is an ik_x multiplier

  !> a(mphiw3,mrad) In the spectral method 
  !> 1st usage: a = grad_y_k <phi_k> = zeta gradient of the gyroaverage of
  !>                potential in k space.
  !> 2nd usage: a = grad_y_k <A||_k>, electromagnetic terms
  !> 3rd usage: a = grad_y_k <B||_k>, compressional terms
  !> 4th usage: a = grad_x f_k = radial gradient of the distribution function
  !>                in k space.
  !> a(mphiw3,-1:nx+2) In the non_spectral method 
  !> 1st usage  a = phi  
  complex, save, allocatable :: a(:,:,:)
  complex, save, allocatable :: aa(:,:,:,:)

  !> b(mphiw3,mrad) In the spectral method 
  !> 1st usage: b = grad_x_k <phi_k> = radial gradient of the gyroaverage of
  !>                potential in k space.
  !> 2nd usage: b = grad_x_k <A||_k>, electromagnetic terms
  !> 3rd usage: a = grad_y_k <B||_k>, compressional terms
  !> 4th usage: b = grad_y f_k = zeta gradient of the distribution function in
  !>                k space.
  !> b(mphiw3,-1:nx+2) in the non spectral method 
  !> 1st usage  b = fdisi 
  complex, save, allocatable :: b(:,:,:)
  complex, save, allocatable :: bb(:,:,:,:) 

  !complex, allocatable :: c(:,:)  !< Another fft dummy array
  complex, save, allocatable :: c(:,:,:)

  ! arrays for precalculated quantities
  complex, save, allocatable :: ci1kxrh(:)
  complex, save, allocatable :: ci1krho(:)
  complex, save, allocatable :: a_phi(:,:,:,:,:)
  complex, save, allocatable :: b_phi(:,:,:,:,:)
  complex, save, allocatable :: a_apar(:,:,:,:,:)
  complex, save, allocatable :: b_apar(:,:,:,:,:)
  complex, save, allocatable :: a_bpar(:,:,:,:,:)
  complex, save, allocatable :: b_bpar(:,:,:,:,:)

  !> In the spectral method ar = grad_y <phi> = zeta gradient of the 
  !> gyroaverage of potential in real space: ar(mphi,mrad)
  !> For the nonspectral method it is the potential ar(mphi,-1:nx+2) 
  !> or aar() 
  real, save, allocatable :: ar(:,:,:)
  real, save, allocatable :: aar(:,:,:,:)

  !> In the spectral method br = grad_x <phi> = radial gradient of the 
  !> gyroaverage of potential in real space. br(mphi,mrad)
  !> For the nonspectral method it is the distribution br(mphi,-1:nx+2)
  real, save, allocatable :: br(:,:,:)
  real, save, allocatable :: bbr(:,:,:,:)

  !> In the spectral method cr = grad_x f = radial gradient of the 
  !> gyroaverage of potential in real space. later reused for the rhs
  !> cr(mphi,mrad)
  !> In the nonspectral method it is used for the rhs cr(mphi,-1:nx+2)
  real, save, allocatable :: cr(:,:,:)
  real, save, allocatable :: ccr(:,:,:,:)

  !> dr = grad_y f = zeta gradient of the gyroaverage of potential in real
  !> space. dr(mphi,mrad)
  real, save, allocatable :: dr(:,:,:)
  real, save, allocatable :: ddr(:,:,:,:)

  !> er = grad_zeta ( 2. vthref  <A||> ) = poloidal gradient of the the
  !> parallel vector potential. er (mphi,mrad)
  real, save, allocatable :: er(:,:,:)
  real, save, allocatable :: eer(:,:,:,:)  


  !> fr = grad_psi ( 2. vthref  <A||> ) = radial gradient of the the parallel
  !> vector potential.  fr (mphi,mrad)
  real, save, allocatable :: fr(:,:,:)
  real, save, allocatable :: ffr(:,:,:,:)

  !> gr = grad_zeta ( 2. mu T_R  <B||> / Z ) = poloidal gradient of bpar
  !> potential. er (mphi,mrad)
  real, save, allocatable :: gr(:,:,:)
  real, save, allocatable :: ggr(:,:,:,:)

  !> hr = grad_psi ( 2. mu T_R  <B||> / Z ) = radial gradient of bpar
  !> potential.  fr (mphi,mrad)
!  real, save, allocatable :: fr(:,:)
  real, save, allocatable :: hr(:,:,:)

  !> The number of 'poloidal' (elsewhere called toroidal) points in the FFT
  !> Notationally this should probably be called mpsi!
  !> this is bigger than nmod because of the dealiasing
  integer, save :: mphi, mphiw3

  !> The number of radial points in the FFT
  !> This is bigger than nx because of the dealiasing
  integer, save :: mrad

  !>Distribution function index translation array
  !>To facilitate copying of the distribution function
  ! integer, allocatable, save :: lindx(:)
  integer, save, allocatable :: lindx3(:,:,:,:), linphi(:,:,:,:)

  ! integer array to copy back the values of the arrays
  ! integer, allocatable, save :: lincopy(:)
  integer, save, allocatable :: lincopy3(:,:,:,:)
  
  !> integer array for the application of the boundary conditions 
  integer, save, allocatable :: ibd(:,:,:) 

  !> real array for the application of the boundary conditions 
  complex, save, allocatable :: mbd(:,:,:) 

  ! number of radial ghost points used in this module (zero unless parallel_x)
  integer, save :: xgp
  ! 
  integer, save :: x_boundary_points

  !> Flag to set once nonlinear_init is called
  logical, save  :: nl_initialised = .false.

  !> maximum number of OpenMP threads
  integer, save  :: max_threads
  integer, save  :: num_threads_request

  !> Index array used in the calculation of the parallel derivative of the 
  !> potential in the velocity nonlinearity ipv(4,ns)  
  integer, save, allocatable :: ipv(:,:)
  
  !> Coefficients used in the calculation of the parallel derivative of the 
  !> potential in the velocity nonlinearity dpvnl(4,ns*nx*nmod) 
  complex, save, allocatable :: dpvnl(:,:) 
  
  !> Index array used in the calculation of the derivative of the distribution 
  !> towards the parallel velocity 
  integer, save, allocatable :: ipvnl(:,:)
  
  !> Coefficient array for the calculation of the derivative of the distribution 
  !> towards the parallel velocity cpvnl(4,nsp*nmu*nvpar*ns*nx*nmod) 
  complex, save, allocatable :: cpvnl(:,:)
  
  !> Coefficient array for the calculation of the drift velocity in connection
  !> with the parallel velocity nonlinearity
  real, save, allocatable :: epvnl(:,:)
  
  !> Coefficient array for species dependent factor in the parallel velocity 
  !> nonlinearity 
  real, save, allocatable :: fpvnl(:)
  
  !> logical that determines which of the two routines for non_spectral is 
  !> called (new one is 2x slower, but more flexible for Sung derivatives)
  logical, save :: nspc2 = .false. 
  
  !> Index array for the implementation of the finite rho* parallel 
  !> derivatives 
  integer, save, allocatable :: ips(:,:) 

  !> Coefficient for the implementation of the finite rho* parallel 
  !> derivatives  
  complex, save, allocatable :: cps(:,:)

  
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> a helper function, returns true if all prime factors of number are smaller
!> or equal to max_prime
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function prime_factors_smallereq_than(number, max_prime)
  integer, intent(in) :: number
  integer, intent(in) :: max_prime
  logical :: prime_factors_smallereq_than
  integer :: i,n
  i = 2
  n = number
  do while(.true.)
    if(mod(n, i) == 0) then
      n = n / i
    elseif(i == max_prime) then
      prime_factors_smallereq_than = (n == 1)
      return
    else
      i = i + 1
    end if
  end do

end function prime_factors_smallereq_than

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_extended_firstdim_fft_size(nmod, posspace_size, kgrid_size)
  integer, intent(in) :: nmod
  integer, intent(out) :: posspace_size, kgrid_size
  integer :: i
  
!#define OLD_GRIDSIZE
#ifdef OLD_GRIDSIZE
  ! Calculate the size of the binormal grid on which the FFT will work.
  real :: dum

  if (nmod == 1) then
    dum = 1.E0
  else
    dum = 1.5*real(2*nmod - 2)
  end if


  posspace_size = log(dum)/log(2.E0) + 1.E0
  posspace_size = 2**posspace_size

#else
  ! As the ExB nonlinearity is quadratic (it goes like phi(k)*f(k)) it
  ! will produce higher modes, up to 2*krhomax, which corresponds to
  ! the array index (2*nmod -1) because imod=1 is the zeromode.

  ! In order to prevent aliasing on the first nmod modes, one needs
  ! kgrid_size > (3/2)*(nmod-1)*2
  ! Then the highest modes are folded back onto that part of the
  ! spectrum which is truncated for dealiasing anyway.

  posspace_size = 3*nmod - 2
  ! make it an even number:
  if(mod(posspace_size, 2) /= 0) then
    posspace_size = posspace_size + 1
  end if
  
  do while(.not.prime_factors_smallereq_than(posspace_size, 7))
    posspace_size = posspace_size + 2
  end do

  ! The division by two is rounded down, according to the FFTW3 manual.
  ! Hence, for every kgrid_size there is one even and one odd value
  ! of posspace_size possible.
  ! Scans of nmod showed that the main loop time was increased in cases where
  ! posspace_size = (kgrid_size-1)*2+1 only has small prime factors, but the
  ! corresponding second number (kgrid_size-1)*2 has at least one large
  ! prime factor.

  ! For this reason, posspace_size is chosen to be even.  The FFTW3
  ! manual points out that it is also beneficial for perfomance to
  ! have an even size of the fastest running dimension (= the first
  ! dimension in Fortran).

  ! According to performance measurements, it could be faster to take
  ! a 2**n grid, if it is not very much larger. So, have a look, if there
  ! is one closeby. 
  look_at_close_alternatives: do i = 1,8
     ! (consider only even numbers)
    if(prime_factors_smallereq_than(posspace_size + i*2, 2)) then
      posspace_size = posspace_size + i*2
      exit look_at_close_alternatives;
    end if
  end do look_at_close_alternatives

#endif

  ! -> posspace_size is the number of points of the binormal grid in position space.
  ! Because of the symmetry of the Fourier representation of
  ! a real-valued function, it is sufficient to store only (posspace_size/2 + 1)
  ! values in k-space.
  ! http://www.fftw.org/fftw3_doc/Real_002ddata-DFT-Array-Format.html#Real_002ddata-DFT-Array-Format
  ! Hence
  ! -> kgrid_size is the number of points of the binormal grid in k-space
  kgrid_size = (floor(posspace_size/2.0) + 1)

end subroutine get_extended_firstdim_fft_size

!-----------------------------------------------------------------------------
!>
!-----------------------------------------------------------------------------
subroutine get_extended_seconddim_fft_size(posspace_size, kgrid_size)
  integer, intent(in) :: posspace_size
  integer, intent(out) :: kgrid_size
  real :: dum
  integer :: i
  ! In order to prevent aliasing, one needs
  !   kgrid_size > (3/2)*posspace_size .
!#define OLD_RAD_GRIDSIZE
#ifdef OLD_RAD_GRIDSIZE
  ! THE OLD WAY TO GET A kgrid_size VALUE:
  dum  = 1.5*real(posspace_size+1)
  kgrid_size = log(dum)/log(2.) + 1.
  kgrid_size = 2**kgrid_size

#else
#ifdef FFT_FFTW3
  ! Get an appropriate fft grid size:
  ! find a value kgrid_size which is > (3/2)*posspace_size and has small prime factors

  ! (This can be used to get efficient non power of 2 sizes; it is not clear
  !  at what point rounding to the next power of 2 would be most efficient
  !  e.g. 63 and 64 are both supposed to be efficient.)
  ! The FFTW3 library is by default
  ! compiled to be most efficient for prime factors <= 7

  ! Get an appropriate fft grid size using a small function.
  dum  = ceiling(1.5*real(posspace_size+1)) + 1
  do while(.not.prime_factors_smallereq_than(int(dum), 7))
    dum = dum + 1
  end do
#endif

  ! According to performance measurements, it could be faster to take
  ! a 2**n grid, if it is not much larger. So, have a look, if there
  ! is one closeby. 
  look_at_close_alternatives: do i = 1,8
    if(prime_factors_smallereq_than(int(dum) + i, 2)) then
      dum = dum + i
      exit look_at_close_alternatives;
    end if
  end do look_at_close_alternatives
  kgrid_size = int(dum)
#endif


end subroutine get_extended_seconddim_fft_size

!-----------------------------------------------------------------------------
!> This routine allocates the help arrays for the FFTs
!-----------------------------------------------------------------------------
subroutine nonlinear_allocate()

  use control,        only : nlapar, nlbpar, spectral_radius
  use grid,           only : nmod, nx, ns, nvpar, nmu, nsp
  use general,        only : gkw_abort
  use ompinterface,   only : ompget_max_threads
  use dist,           only : ghost_points_xf, stencil_side
  use rho_par_switch, only : lnonlinear_rhostar
  use control,        only : lpar_vel_nl
  use global, only : id_x

  integer :: ierr

  ! The Sung terms and parallel velocity nonlinearity are only implemented 
  ! in the nspc2 version 
  if (lnonlinear_rhostar) nspc2 = .true.
  if (lpar_vel_nl)        nspc2 = .true. 
  
  ! Get number of openmp threads
  max_threads = ompget_max_threads()

  ! The logical size of position space grids is mphi and mrad (or at least
  ! almost).  These numbers thus determine the resolution and affect also
  ! the nonlinear timestep estimate.
  ! First determine the size of the extended binormal grid. For the binormal
  ! grid, the physical and the binormal grid size differ.
  call get_extended_firstdim_fft_size(nmod, mphi, mphiw3)
  
  ! Now the size of the extended radial grid has to be determined.
  ! NOTE: The symmetry of the Fourier representation can only be used
  ! *in the 1st dimension* to store only half of the spectrum. Here
  ! this is the binormal direction. The 2nd dimension (here 'radial')
  ! cannot profit from this symmetry (cf. FFTW3 docs).

  call get_extended_seconddim_fft_size(nx, mrad)
  ! write (*,*) "mphi:", mphi, "mphiw3:", mphiw3, "mrad:", mrad

  ! set the status integer
  ierr = 0

  xgp = ghost_points_xf
  x_boundary_points = stencil_side(id_x)

  
  ! always used, by diagnostics
  allocate(jind(nx), stat = ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate jind in non_linear_terms')
  end if
 
  if (nspc2 .and. (.not.spectral_radius)) then 
    allocate(lindx3(nvpar*nmod*(nx+2*xgp)*ns,1,nmu,nsp),stat = ierr)
  else  
    allocate(lindx3(nvpar*nmod*(nx+2*xgp),ns,nmu,nsp),stat=ierr)
  endif
  if (ierr /= 0) then
    call gkw_abort('Could not allocate lindx in non_linear_terms')
  end if

  if ((.not.spectral_radius).and.nspc2) then 
  allocate(lincopy3(nvpar*nmod*nx*ns,1,nmu,nsp),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate lincopy in non_linear_terms')
  end if
  else  
  allocate(lincopy3(nvpar*nmod*nx,ns,nmu,nsp),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate lincopy in non_linear_terms')
  end if
  endif

  allocate(i3loop(3,nsp*nmu*ns),stat=ierr)
  if (ierr /= 0) then
    call gkw_abort('Could not allocate i3loop in non_linear_terms')
  end if
  

  spectral: if (spectral_radius) then 

    allocate(ci1krho(nmod),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate ci1krho in non_linear_terms')
    end if

    allocate(ci1kxrh(nx),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate ci1krho in non_linear_terms')
    end if

    allocate(a(mphiw3,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate a in non_linear_terms')
    end if

    allocate(b(mphiw3,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate b in non_linear_terms')
    end if

    allocate(ar(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate ar in non_linear_terms')
    end if

    allocate(br(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate br in non_linear_terms')
    end if

    allocate(cr(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate cr in non_linear_terms')
    end if

    allocate(dr(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate dr in non_linear_terms')
    end if

    allocate(er(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate er in non_linear_terms')
    end if

    allocate(fr(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate fr in non_linear_terms')
    end if

    allocate(gr(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate gr in non_linear_terms')
    end if

    allocate(hr(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate hr in non_linear_terms')
    end if

    allocate(c(mphi,mrad,max_threads), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate c in non_linear_terms')
    end if

    allocate(jinv(mrad), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate jinv in non_linear_terms')
    end if

    allocate(a_phi(nmod,nx,ns,nmu,nsp),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate a_phi in non_linear_terms')
    end if

    allocate(b_phi(nmod,nx,ns,nmu,nsp),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate b_phi in non_linear_terms')
    end if

    if (nlapar) then
      allocate(a_apar(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate a_apar in non_linear_terms')
      end if

      allocate(b_apar(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate b_apar in non_linear_terms')
      end if
    end if

    if (nlbpar) then
      allocate(a_bpar(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate a_bpar in non_linear_terms')
      end if

      allocate(b_bpar(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate b_bpar in non_linear_terms')
      end if
    end if

  else 
   
    nonspec_method: if (nspc2) then 

      allocate(aa(mphiw3,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate a in non_linear_terms')

      allocate(bb(mphiw3,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate b in non_linear_terms')

      allocate(aar(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate ar in non_linear_terms')

      allocate(bbr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate br in non_linear_terms')

      allocate(ccr(mphi,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate cr in non_linear_terms')

      allocate(ddr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate dr in non_linear_terms')

      allocate(ibd(-1:2,nmod,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate ibd in non_linear terms')

      allocate(mbd(-1:2,nmod,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate mbd in non_linear terms')

      if (nlapar) then

        allocate(eer(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
        if (ierr /= 0) call gkw_abort('Could not allocate er in non_linear_terms')

        allocate(ffr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads), stat = ierr)
        if (ierr /= 0) call gkw_abort('Could not allocate fr in non_linear_terms')

      end if

      if (nlapar) then
        allocate(linphi(2*nmod*(nx+2*xgp)*ns,1,nmu,nsp), stat = ierr)
      else
        allocate(linphi(nmod*(nx+2*xgp)*ns,1,nmu,nsp), stat = ierr)
      end if
      if (ierr /= 0) call gkw_abort('Could not allocate linphi in non_linear_terms')

      if (lnonlinear_rhostar) then 

        allocate(ips(ns,-2:2), stat = ierr) 
        if (ierr /= 0) call gkw_abort('Could not allocate ips in non_linear_terms')
        allocate(cps(ns,-2:2), stat = ierr) 
        if (ierr /= 0) call gkw_abort('Could not allocate cps in non_linear_terms')
        allocate(ggr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,ns,max_threads),stat = ierr) 
        if (ierr /= 0) call gkw_abort('Could not allocate ggr in non_linear_terms')

      endif

      if (lpar_vel_nl) then 

        allocate(ipv(4,ns), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate ipv in non_linear_terms')
        allocate(ipvnl(4,nsp*nmu*nvpar*ns*nx*nmod), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate ipvnl in non_linear_terms')
        allocate(cpvnl(4,nsp*nmu*nvpar*ns*nx*nmod), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate cpvnl in non_linear_terms')
        allocate(dpvnl(4,ns*nmod*nx), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate dpvnl in non_linear_terms')
        allocate(epvnl(2,nsp*nmu*nvpar*ns*nx*nmod), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate epvnl in non_linear_terms')
        allocate(fpvnl(nsp), stat = ierr) 
        if (ierr/=0) call gkw_abort('Could not allocate fpvnl in non_linear_terms') 

      endif

    else 

      allocate(a(mphiw3,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate a in non_linear_terms')

      allocate(b(mphiw3,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate b in non_linear_terms')

      allocate(ar(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate ar in non_linear_terms')

      allocate(br(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate br in non_linear_terms')

      allocate(cr(mphi,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate cr in non_linear_terms')

      allocate(dr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate dr in non_linear_terms')

      allocate(ibd(-1:2,nmod,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate ibd in non_linear terms')

      allocate(mbd(-1:2,nmod,ns), stat = ierr) 
      if (ierr /= 0) call gkw_abort('Could not allocate mbd in non_linear terms')

      if (nlapar) then

        allocate(er(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
        if (ierr /= 0) call gkw_abort('Could not allocate er in non_linear_terms')

        allocate(fr(-1:mphi+2,1-x_boundary_points:nx+x_boundary_points,max_threads), stat = ierr)
        if (ierr /= 0) call gkw_abort('Could not allocate fr in non_linear_terms')

      end if !nlapar

      if (nlapar) then
        allocate(linphi(2*nmod*(nx+2*xgp),ns,nmu,nsp), stat = ierr)
      else
        allocate(linphi(nmod*(nx+2*xgp),ns,nmu,nsp), stat = ierr)
      end if
      if (ierr /= 0) call gkw_abort('Could not allocate linphi in non_linear_terms')

    endif nonspec_method

  endif spectral

end subroutine nonlinear_allocate

!-----------------------------------------------------------------------------
!> This routine performs initialisations of quantites
!> Note it may be called even if add_non_linear terms is not,
!> As the parameters it calculates are in used for mode_box,
!> 2D diagnostics and some of the experiemental shear routines.
!>
!> Sets up indexing arrays, calculates FFT size.
!-----------------------------------------------------------------------------
subroutine nonlinear_init()
  use control,        only : nlapar, nlbpar, non_linear, spectral_radius 
  use control,        only : shift_metric, vp_trap, lpar_vel_nl
  use constants,      only : ci1
  use mode,           only : krho, kxrh, ixzero, parallel_phase_shift, lyinv
  use index_function, only : indx
  use dist,           only : ifdis, iphi_ga, iapar_ga
  use geom,           only : ffun, dfun, efun, hfun, bn, sgr_dist
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp, parallel_s, lx 
  use general,        only : gkw_abort, gkw_warn
  use fft,            only : four2D_real, four1D_real, FFT_FORWARD, FFT_INVERSE
  use functions,      only : besselj0_gkw, mod_besselj1_gkw
  use components,     only : vthrat, tmp, signz, tgrid, mas, rhostar
  use components,     only : veta_prime
  use velocitygrid,   only : mugr, vpgr
  use global,         only : compiled_with_openmp, int2char
  use structures,     only : matrix_element
  use matdat,         only : connect_rad, set_indx
  use rho_par_switch, only : lnonlinear_rhostar
  use mpiinterface,   only : root_processor
  use rotation,       only : vcor
  use ompinterface,   only : ompget_max_threads
  ! Just loop indices
  integer :: imod, ix, i, j, k, jv, kt, is, ipar, idx, idx2, idum, inp
  real    :: b0, b1_mod, ffac = 0.0
  complex :: cdum = (0.0, 0.0)
  logical :: ingrid 
  type (matrix_element) :: E 

  ! This initialization routine should only be called once
  if (nl_initialised) call gkw_abort('nonlinear_init called twice')

  E%term = ""

  if (root_processor .and. compiled_with_openmp .and. non_linear) then
    idum = nsp*ns*nmu
    if ((.not. spectral_radius) .and. nspc2) idum = nsp*nmu
    write(*,*)
    ! it depends on some env variables how many threads actually get, but no
    ! more than this:
    write(*,*) 'OpenMP threads requested for NL terms computation:', max_threads
    write(*,*) 'NL outer loop iterations (given by local &
       & gridsize nsp*ns*nmu):', idum
    write(*,*)
  end if
  
  ! Check openmp setup
  if (compiled_with_openmp .and. mod(nsp*ns*nmu,max_threads) /= 0 .and. &
     & non_linear) then   
    ! search a better number of threads
#ifdef _OPENMP
    recommend_num_threads: do i = ompget_max_threads(), 1, -1
      if (mod(nsp*ns*nmu, i) == 0) then
        call gkw_warn('Requested number of threads should factor into NL&
           & outer loop iterations, for most effective use of the ''static''&
           & schedule. Hence e.g. '// &
           & int2char(i)//' threads would be better.')
        num_threads_request = max_threads
        exit recommend_num_threads
      end if
    end do recommend_num_threads
#endif
  else
    ! if the implementation is not capable of supporting the requested
    ! number of threads, the behaviour is implementation defined.
    num_threads_request = max_threads
  end if
  
  if (compiled_with_openmp .and. max_threads > 1 .and. .not. spectral_radius) then
    call gkw_warn('Open MP not implemented in non-spectral fields solve')
  end if

  if (spectral_radius) then
    ! The array over which the radial wavevectors are stored differs
    ! in the code and in the fft.  The jind array translates.
    ! For fast evaluation the indices are calculated
    jind = 0
    do ix = ixzero, nx
      jind(ix) = ix - ixzero + 1
    end do
    do ix = ixzero-1, 1, -1
      jind(ix) = mrad + ix - ixzero + 1
    end do
    
    ! The jinv array provides the inverse translation
    do i = 1, mrad
      jinv(i) = 0
      do j = 1, nx
        if (jind(j) == i) jinv(i) = j
      end do
    end do
  else
     ! for the benefit of simplicity 
    do ix = 1, nx
      jind(ix) = ix
    end do     
  endif

  ! write (*,*) "ixzero:", ixzero
  ! write (*,*) "jind:", jind
  ! write (*,*) "jind_flexible:", (/ (jind_flexible(nx, mrad, i), i = 1,nx)/)
  ! write (*,*) "jinv:", jinv
  ! write (*,*) "jinv_flexible:", (/ (jinv_flexible(nx, mrad, i), i = 1,mrad)/)

  ! construct the array for quick look up of the distribution
  lindx3 = 0

  if (nspc2.and.(.not.spectral_radius)) then 

    do is = 1, nsp; do jv = 1, nmu 
      idx2 = 1 
      do kt = 1, nvpar; do ipar = 1, ns; do ix = 1-xgp, nx+xgp; do imod = 1, nmod
        lindx3(idx2,1,jv,is) = indx(ifdis,imod,ix,ipar,jv,kt,is)
        idx2 = idx2 + 1 
      end do; end do; end do; end do; 
    end do; end do 
 
  else 

    do is = 1, nsp ; do jv = 1, nmu ; do ipar = 1, ns
     idx2 = 1
      do kt = 1, nvpar
        do ix = 1-xgp, nx+xgp
          do imod = 1, nmod
            lindx3(idx2,ipar,jv,is) = indx(ifdis,imod,ix,ipar,jv,kt,is)
            idx2 = idx2 + 1
          end do
        end do
      end do
    end do ; end do ; end do

  endif 

  ! Construct indexing array for Collapsing OpenMP (!!$OMP) nested loops
  idx = 1
  do is = 1, nsp         ! Loop over species
    do jv = 1, nmu       ! Loop over magnetic moment
      do ipar = 1, ns         ! Loop along field line points
        i3loop(1,idx) = is
        i3loop(2,idx) = jv
        i3loop(3,idx) = ipar
        idx = idx + 1
      end do
    end do
  end do

  if (.not. spectral_radius) then
  
    ! construct the array for quick look up of the ga fields
    linphi = 0
 
    if (nspc2) then 

      ! also recalculate i3loop 
      idx = 1 
      do is = 1, nsp; do jv = 1, nmu; 
        i3loop(1,idx) = is 
        i3loop(2,idx)  = jv 
        idx = idx + 1 
      end do; end do; 

      do is = 1, nsp ; do jv = 1, nmu 
        idx2=1
        ! do apar field first for copy convenience
        if (nlapar) then
          do ipar = 1, ns; do ix = 1-xgp, nx+xgp; do imod = 1, nmod 
             linphi(idx2,1,jv,is) = indx(iapar_ga,imod,ix,ipar,jv,is)
             idx2 = idx2 + 1
          end do; end do; end do 
        end if 
        do ipar = 1, ns; do ix = 1-xgp, nx+xgp; do imod = 1, nmod 
          linphi(idx2,1,jv,is) = indx(iphi_ga,imod,ix,ipar,jv,is)
          idx2 = idx2 + 1
        end do; end do; end do 
      end do ; end do 
 
    else 
   
      do is = 1, nsp ; do jv = 1, nmu ; do ipar = 1, ns
        idx2=1
        ! do apar field first for copy convenience
        if (nlapar) then
          do ix = 1-xgp, nx+xgp; do imod = 1, nmod
             linphi(idx2,ipar,jv,is) = indx(iapar_ga,imod,ix,ipar,jv,is)
             idx2 = idx2 + 1
          end do; end do
        end if

        do ix = 1-xgp, nx+xgp; do imod = 1, nmod
          linphi(idx2,ipar,jv,is) = indx(iphi_ga,imod,ix,ipar,jv,is)
          idx2 = idx2 + 1
        end do; end do
      end do ; end do ; end do

    endif 
    
    ! construct the array for copy back
    ! array to copy back the nonlinear terms in the right hand side
    lincopy3 = 0

    if (nspc2) then 

      do is = 1, nsp ; do jv = 1, nmu 
        idx2=1
        do kt = 1, nvpar; do j = 1, nx; do i = 1, nmod; do ipar = 1, ns
          lincopy3(idx2,1,jv,is) = indx(ifdis,i,j,ipar,jv,kt,is)
          idx2 = idx2 + 1
        end do ; end do; end do; end do 
      end do ; end do 

    else 

      do is = 1, nsp ; do jv = 1, nmu ; do ipar = 1, ns
        idx2=1
        do kt = 1, nvpar; do j = 1, nx; do i = 1, nmod
          lincopy3(idx2,ipar,jv,is) = indx(ifdis,i,j,ipar,jv,kt,is)
          idx2 = idx2 + 1
        end do ; end do; end do
      end do ; end do ; end do

    endif

    ! Call the FFTs once to force FFTW plan setup before OMP section
    ! The planning routines are not threadsafe !
    if (nspc2) then 

      aa(:,1,1,1) = (0.0,0.0)
      aar(:,1,1,1) = 0.0
      call four1D_real(aar(1:mphi,1,1,1),aa(:,1,1,1),FFT_INVERSE)
      call four1D_real(aar(1:mphi,1,1,1),aa(:,1,1,1),FFT_FORWARD)

    else 

      a(:,1,1) = (0.0,0.0)
      ar(:,1,1) = 0.0
      call four1D_real(ar(1:mphi,1,1),a(:,1,1),FFT_INVERSE)
      call four1D_real(ar(1:mphi,1,1),a(:,1,1),FFT_FORWARD)

    endif

    ! Set up the boundary conditions
    do imod = 1, nmod; do i = 1, ns

      ! set a dummy index to 1 and set the index
      idum = 1 
      call set_indx(E,imod,idum,i,idum,idum,idum)
      E%val   = 1.
      E%ixloc = -1
      call connect_rad(E,ingrid)
      ibd(-1,imod,i) = E%ixloc
      mbd(-1,imod,i) = E%val

      ! set a dummy index to 1 and set the index
      idum = 1 
      call set_indx(E,imod,idum,i,idum,idum,idum)
      E%val   = 1.
      E%ixloc = 0
      call connect_rad(E,ingrid)
      ibd(0,imod,i) = E%ixloc
      mbd(0,imod,i) = E%val

      ! set a dummy index to 1 and set the index
      idum = 1
      call set_indx(E,imod,idum,i,idum,idum,idum)
      E%val   = 1.
      E%ixloc = nx+1
      call connect_rad(E,ingrid)
      ibd(1,imod,i) = E%ixloc
      mbd(1,imod,i) = E%val

      ! set a dummy index to 1 and set the index
      idum = 1
      call set_indx(E,imod,idum,i,idum,idum,idum)
      E%val   = 1.
      E%ixloc = nx+2
      call connect_rad(E,ingrid)
      ibd(2,imod,i) = E%ixloc
      mbd(2,imod,i) = E%val

    end do; end do 

    if (nspc2.and.lnonlinear_rhostar) then

      ! the implementation does not work with parallel_s or
      ! with shifted metric
      if (parallel_s .or. shift_metric) &
        & call gkw_abort('NL Sung does not work with parallel_s or shift_metric')

      do i = 1, ns; do j = -2,2
        ips(i,j) = i + j
        if (ips(i,j).lt.1)  ips(i,j) = ns + ips(i,j)
        if (ips(i,j).gt.ns) ips(i,j) = ips(i,j) - ns
        cps(i,j) = parallel_phase_shift(1,1,1,i+j)
      end do; end do
      do i = 1, ns
        cps(i,-2) =  1.0E0*cps(i,-2)
        cps(i,-1) = -8.0E0*cps(i,-1)
        cps(i, 0) =  0.0E0*cps(i, 0)
        cps(i, 1) = +8.0E0*cps(i, 1)
        cps(i, 2) = -1.0E0*cps(i, 2)
      end do

    endif

    if (nspc2.and.lpar_vel_nl) then 

      ! Warning loop over the fields might cause problems. Run only with 
      ! electro-static at the moment 

      ! set the counter 
      inp = 0 
      
      ! First the arrays that are used for the derivative towards the parallel 
      ! velocity (The boundary conditions are simplified here, and perhaps 
      ! need revisiting) 

      do is = 1, nsp; 

        ! The normalizing factor. This is only valid for vp_trap = 0. 
        ! test on the vp_trap should be moved elsewhere 
        if (vp_trap == 0) then  
          ffac = 1.E0 / (12.E0 * (vpgr(1,1,2,is)-vpgr(1,1,1,is)))
        else 
          call gkw_abort('Parallel velocity nonlinearity can only be run with & 
                       & a uniform parallel velocity grid (vp_trap = 0)')
        endif       
      
        ! a common normalization factor is added to the parallel velocity 
        ! derivative 
        do j = 1, nmu; do k = 1, nvpar; do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
          inp = inp + 1 
          kt = k - 2 
          if (kt.ge.1) then 
            ipvnl(1,inp) = indx(ifdis,imod,ix,i,j,kt,is)
            cpvnl(1,inp) = +1.E0 * rhostar * ffac 
          else 
            ipvnl(1,inp) = indx(ifdis,imod,ix,i,j,k,is)
            cpvnl(1,inp) = +0.E0 * rhostar * ffac 
          endif             
          kt = k - 1 
          if (kt.ge.1) then 
            ipvnl(2,inp) = indx(ifdis,imod,ix,i,j,kt,is)
            cpvnl(2,inp) = -8.E0 * rhostar * ffac
          else 
            ipvnl(2,inp) = indx(ifdis,imod,ix,i,j,k,is)
            cpvnl(2,inp) = +0.E0 * rhostar * ffac
          endif             
          kt = k + 1 
          if (kt.le.nvpar) then 
            ipvnl(3,inp) = indx(ifdis,imod,ix,i,j,kt,is)
            cpvnl(3,inp) = +8.E0 * rhostar * ffac
          else 
            ipvnl(3,inp) = indx(ifdis,imod,ix,i,j,k,is)
            cpvnl(3,inp) = +0.E0 * rhostar * ffac
          endif             
          kt = k + 2 
          if (kt.le.nvpar) then 
            ipvnl(4,inp) = indx(ifdis,imod,ix,i,j,kt,is)
            cpvnl(4,inp) = -1.E0 * rhostar * ffac
          else 
            ipvnl(4,inp) = indx(ifdis,imod,ix,i,j,k,is)
            cpvnl(4,inp) = +0.E0 * rhostar * ffac
          endif             
        end do; end do; end do; end do; end do
      end do 

      ! ffac is a common factor that appears in the parallel derivative contribution 
      ! it is Z / (2.0 sqrt(T_G m)) 
      do is = 1, nsp 
        fpvnl(is) = signz(is) / (2.0E0 * sqrt(tgrid(is)*mas(is))) 
      end do 
      
      ! set the counter 
      inp = 0 

      ! arrays dpvnl and ipv are used in the calculation of the parallel derivative of 
      ! the potential 
      do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
        inp = inp + 1 
        do j = 1,4
          select case (j) 
          case (1) ; ipv(j,i) = i - 2 ; cdum = +1.0E0
          case (2) ; ipv(j,i) = i - 1 ; cdum = -8.0E0
          case (3) ; ipv(j,i) = i + 1 ; cdum = +8.0E0
          case (4) ; ipv(j,i) = i + 2 ; cdum = -1.0E0
          end select 
          ! This coefficient contains several factors 
          ! cdum * ffun / 12.0 * sgr_dist is the factor that follows from the parallel 
          !                               derivative 
          ! parallel_phase_shift is due to the boundary condition at the end of the field 
          !                      line 
          ! mphi is the normalizing coefficient for the inverse Fourier transform (Equal 
          !      to the number of points in the FFT)
          dpvnl(j,inp) = cdum*ffun(ix,i)*parallel_phase_shift(imod,ix,i,ipv(j,i))  &
                       & / (12.E0 * sgr_dist * mphi ) 
          ! mapping back the points over the boundary. 
          if (ipv(j,i).lt.1)  ipv(j,i) = ns + ipv(j,i) 
          if (ipv(j,i).gt.ns) ipv(j,i) = ipv(j,i) - ns 
        end do 
      end do; end do; end do  

      ! set the counter 
      inp = 0 
      
      ! The drift components such that pvnl = epvnl * grad phi_N d f / d v_|| / rhostar
      do is = 1, nsp; do j = 1, nmu; do k = 1, nvpar; 
        do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
          inp = inp + 1 
          epvnl(1,inp) = 0.5E0*(vpgr(i,j,k,is)*dfun(ix,i,1)                               &
                       &     +  vpgr(i,j,k,is)*veta_prime(ix)*efun(ix,i,1,1)/bn(ix,i)**2  &
                       &     +  2.E0*hfun(ix,i,1)*vcor*sqrt(mas(is)/tgrid(is)))  
          epvnl(2,inp) = 0.5E0*(vpgr(i,j,k,is)*dfun(ix,i,2)                               &
                       &     +  vpgr(i,j,k,is)*veta_prime(ix)*efun(ix,i,1,2)/bn(ix,i)**2  &
                       &     +  2.E0*hfun(ix,i,2)*vcor*sqrt(mas(is)/tgrid(is)))  
        end do; end do; end do; 
      end do; end do; end do 
      
      ! multiply with the factor 1/(12 delta x/xi) An extra factor 1/mphi appears because 
      ! of the fourier transform inverse fourier transform  
      do i = 1, inp   
        epvnl(1,i) = epvnl(1,i) * lyinv / 12.E0 
        epvnl(2,i) = epvnl(2,i) * mrad / (12.E0* lx * mphi) 
      end do 
      
    endif 
    
    nl_initialised = .true.
    return

  end if  !spectral_radius

  ! precalculated quantities: complex krho and kxrh
  ci1krho(:) = ci1*krho(:)
  ci1kxrh(:) = ci1*kxrh(:)

  ! array to copy back the nonlinear terms in the right hand side
  lincopy3 = 0

  do is = 1, nsp ; do jv = 1, nmu ; do ipar = 1, ns

    idx2=1
    do kt = 1, nvpar

      do j = 1, nx-ixzero+1
        do i = 1, nmod
          lincopy3(idx2,ipar,jv,is) = indx(ifdis,i,jinv(j),ipar,jv,kt,is)
          idx2 = idx2 + 1
        end do
      end do

      do j = mrad+2-ixzero, mrad
        do i = 1, nmod
          lincopy3(idx2,ipar,jv,is) = indx(ifdis,i,jinv(j),ipar,jv,kt,is)
          idx2 = idx2 + 1
        end do
      end do

    end do

  end do ; end do ; end do

  ! precalculated quantities involving bessel functions
  do is=1,nsp ; do jv=1,nmu ; do ipar=1,ns ; do ix=1,nx ; do imod=1,nmod

    b0 = besselj0_gkw(imod,ix,ipar,jv,is)
    a_phi(imod,ix,ipar,jv,is) = ci1krho(imod)*b0
    b_phi(imod,ix,ipar,jv,is) = ci1kxrh(ix)*b0

    if (nlapar) then
      a_apar(imod,ix,ipar,jv,is) = 2.*vthrat(is)*ci1krho(imod)*b0
      b_apar(imod,ix,ipar,jv,is) = 2.*vthrat(is)*ci1kxrh(ix)*b0
    end if

    if (nlbpar) then
      b1_mod = mod_besselj1_gkw(imod,ix,ipar,jv,is)/signz(is)
      a_bpar(imod,ix,ipar,jv,is) = 2.*mugr(jv)*tmp(ix,is)*b1_mod*ci1krho(imod)
      b_bpar(imod,ix,ipar,jv,is) = 2.*mugr(jv)*tmp(ix,is)*b1_mod*ci1kxrh(ix)
    end if

  end do       ; end do      ; end do      ; end do     ; end do

  ! Call the FFTs once to force FFTW plan setup before OMP section
  ! The planning routines are not threadsafe !
  ! What happens if this is called with nofft?
  a(:,:,1) = (0.0,0.0)
  dr(:,:,1) = 0.0
  call four2D_real(dr(:,:,1),a(:,:,1),FFT_INVERSE)
  call four2D_real(dr(:,:,1),a(:,:,1),FFT_FORWARD)

  ! set initialised to true
  nl_initialised = .true.

end subroutine nonlinear_init

!------------------------------------------------------------------------------
!> Term III in the manual
!------------------------------------------------------------------------------
subroutine add_non_linear_terms(fdis,rhs)

  use dist,    only : nsolc 
  use control, only : spectral_radius
  use general, only : gkw_warn

  complex, intent(in)    :: fdis(:)
  complex, optional, intent(inout) :: rhs(nsolc)

  if (present(rhs)) then
    ! the nonlinear terms is calculated as part of the right hand side
    ! of the gyrokinetic equation
    if (spectral_radius) then
      call add_non_linear_terms_spectral(fdis,rhs)
    else
      if (nspc2) then
        call add_non_linear_terms_nonspec_2(fdis,rhs) 
      else 
        call add_non_linear_terms_nonspec(fdis,rhs) 
      endif
    endif
  else
    ! the nonlinear terms calculation is invoked by diagnostics
    if (spectral_radius) then
      call add_non_linear_terms_spectral(fdis)
    else
      if(nspc2) then
        ! never happens if diagnostics do their checks properly
        call gkw_warn("diagnostics called add_non_linear_terms ver2 in &
           & nonspectral. (Not implemented)")
      else
        call add_non_linear_terms_nonspec(fdis)
      end if
    end if
  end if

end subroutine add_non_linear_terms

!-----------------------------------------------------------------------------
!> The routine for the spectral method 
!-----------------------------------------------------------------------------
subroutine add_non_linear_terms_spectral(fdis,rhs)

  use dist,           only : nsolc, get_phi, phi, get_apar, apar, get_bpar, bpar
  use mode,           only : lxinv, lyinv, ixzero
  use mode,           only : erase_any_transfer_to, erase_any_transfer_from
  use geom,           only : efun
  use control,        only : dtim, dtim_est, dtim_est_save
  use control,        only : non_linear, nlapar, nl_dtim_est, nlbpar
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp, vpmax
  use grid,           only : gs
  use velocitygrid,   only : vpgr
  use fft,            only : four2D_real, FFT_INVERSE, FFT_FORWARD
  use general,        only : gkw_abort
  use global,         only : r_tiny, r_huge
  use constants,      only : c1
  use rotation,       only : shear_real, grad_pot
  use ompinterface,   only : ompget_thread_num, update_thread_statistics
  use mpiinterface,   only : mpiallreduce_min, number_of_processors
  use perform,        only : perfon, perfoff, perf_measure
  use matdat,         only : get_f_from_g
  use functions,      only : besselj0_gkw

  complex, intent(in)    :: fdis(nsolc)
  complex, optional, intent(inout) :: rhs(nsolc)

  complex :: dtim_cmplx, mphimrad1, cdum
  integer :: idx, idxcopy, idx2, idxcopy2
  integer :: imod, ix, i, j, jv, kt, is, ipar, iloop, th
  real    :: dum, maxvalue, dtim_est_dum, maxvalapar, maxvalbpar
  real    :: dtim_est_apar

  logical :: rhs_perf_measure

  rhs_perf_measure = present(rhs) .and. perf_measure

  ! use complex dtim
  dtim_cmplx = c1*dtim

  ! use mphimrad1i =  c1 / (mrad*mphi)
  mphimrad1 = c1 / (mrad*mphi)

  ! Abort if called incorrectly
  if (.not. (non_linear .or. shear_real)) then
    call gkw_abort('Invalid call to add_non_linear_terms')
  end if

  if (.not. nl_initialised) then
    call gkw_abort('add_non_linear terms: cannot call before nonlinear_init')
  end if

  ! Performance monitoring by timing sections of code
  if (rhs_perf_measure) call perfon ('nlmain',1)
  if (rhs_perf_measure) call perfon ('obtain fields',1)

  ! set the value for the estimate of the timestep to zero (The restriction
  ! due to the electromagnetic scheme are measured separately with maxvalapar)
  maxvalue   = 0.
  maxvalapar = 0.
  maxvalbpar = 0.
  dtim_est_apar = dtim + 1.0

  ! update the fields in the phi, apar and bpar arrays from dist
  if (non_linear) then
    call get_phi(fdis,phi)
    if (nlapar) call get_apar(fdis,apar)
    if (nlbpar) call get_bpar(fdis,bpar)
  end if
  if (rhs_perf_measure) call perfoff(1)

  ! initialize the indices of the help arrays
  idx = 1
  idxcopy = 1

  !!FJC 25-03-2010 OMP TO DO:
  !0) 19-10-2010: Test the suggested index function for parallel MPI/OMP runs.
  !1) Reorder code so a,c,e,g and b,d,f,h loops and FFTS are done in separate
  !   blocks
  !2) Is this 1) slower for single thread ?
  !3) Use $OMP SECTIONS to get effcient scaling up to 4 threads at maximum of
  !   MPI (2 in outer loop)

  ! loop_ns: do ipar = 1, ns           ! Loop along field line points
  !   loop_nsp: do is = 1, nsp         ! Loop over species
  !     loop_nmu: do jv = 1, nmu       ! Loop over magnetic moment

  ! Manually collapse the three nested loops above
!!$OMP do collapse(3) Not supported until OpenMP 3.0

  !$OMP parallel default(shared) num_threads(num_threads_request) &
  !$OMP private(ipar,iloop,is,jv,ix,imod,i,j,kt,dum,th,idx2, &
  !$OMP idxcopy2,cdum)
  !$OMP master
  call update_thread_statistics()
  !$OMP end master

  !$OMP do &
  !$OMP reduction(max:maxvalue,maxvalapar,maxvalbpar)
  loop3: do iloop = 1, ns*nsp*nmu
    is = i3loop(1,iloop)
    jv = i3loop(2,iloop)
    ipar = i3loop(3,iloop)

  ! !$omp do &
  ! !$omp collapse(3) reduction(max:maxvalue,maxvalapar,maxvalbpar)
  ! loop_nsp: do is = 1, nsp
  !   loop_nmu: do jv = 1, nmu
  !     loop_ns: do ipar = 1, ns

    ! Always 1 for non openmp runs
    th = ompget_thread_num() + 1

      ! fill the array for the potential
      get_real_grad_phi: if(non_linear) then

        if (rhs_perf_measure) call perfon('k-space multiplication 1',1)

        ! Gyroaveraged potential in k space obtained with Bessel function J_0
        ! <phi_k> = J_0(k_perp rho) phi_k
        ! Does not depend on parallel velocity, calculated outside loop

        ! a=grad_y_k <phi_k> = zeta gradient of the gyroaverage of potential
        !                      = i J_0() k_zeta phi_k
        a(:,:,th) = (0.,0.)

        ! b=grad_x_k <phi_k> = radial gradient of the gyroaverage of potential
        !                     = i J_0() k_psi phi_k
        b(:,:,th) = (0.,0.)

        loop_nx1: do ix = 1, nx
          loop_nmod1: do imod = 1, nmod
            a(imod,jind(ix),th) = a_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
            b(imod,jind(ix),th) = b_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
          end do loop_nmod1
        end do  loop_nx1

        if (rhs_perf_measure) call perfoff(1)
        if (rhs_perf_measure) call perfon('fft 1',1)

        ! Not openmp threadsafe ?
        if(erase_any_transfer_from) then
          call erase_transfer_from_modes(a,th)
          call erase_transfer_from_modes(b,th)
        end if

        ! Inverse fourier transform of potential k space to real space
        ! ar = grad_zeta <phi> =
        ! poloidal gradient of the gyroaverage of potential in real space
        ! br = grad_psi <phi> =
        ! radial gradient of the gyroaverage of potential in real space

        ! Could subdivide this one transform per thread, if doing the nested
        ! parallism vpar loop below.
        call four2D_real(ar(:,:,th),a(:,:,th),FFT_INVERSE)
        call four2D_real(br(:,:,th),b(:,:,th),FFT_INVERSE)

        ! The electro-magnetic corrections if needed
        if (nlapar) then

          ! initialization
          a(:,:,th) = (0.,0.)
          b(:,:,th) = (0.,0.)

          do ix = 1, nx         ! Loop over radial modes
            do imod = 1, nmod     ! Loop over poloidal modes
              a(imod,jind(ix),th) =                                            &
                 &    a_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
              b(imod,jind(ix),th) =                                            &
                 &    b_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
            end do
          end do

          if(erase_any_transfer_from) then
            call erase_transfer_from_modes(a,th)
            call erase_transfer_from_modes(b,th)
          end if


          ! inverse fourier transform (to real space)
          call four2D_real(er(:,:,th),a(:,:,th),FFT_INVERSE)
          call four2D_real(fr(:,:,th),b(:,:,th),FFT_INVERSE)

        end if

        if (nlbpar) then

          ! initialization
          a(:,:,th) = (0.,0.)
          b(:,:,th) = (0.,0.)

          do ix = 1, nx         ! Loop over radial modes
            do imod = 1, nmod     ! Loop over poloidal modes
              a(imod,jind(ix),th) =                                            &
                 &    a_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
              b(imod,jind(ix),th) =                                            &
                 &    b_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
            end do
          end do

          if(erase_any_transfer_from) then
            call erase_transfer_from_modes(a,th)
            call erase_transfer_from_modes(b,th)
          end if


          ! inverse fourier transform (to real space)
          call four2D_real(gr(:,:,th),a(:,:,th),FFT_INVERSE)
          call four2D_real(hr(:,:,th),b(:,:,th),FFT_INVERSE)

        end if

        if (rhs_perf_measure) call perfoff(1)
        if (rhs_perf_measure) call perfon('estimate max velocity',1)

        ! Estimate maximum velocity for timestep estimator
        maxvalue=max(maxval(abs(ar(:,:,th)))*mrad * lxinv ,maxvalue)
        maxvalue=max(maxval(abs(br(:,:,th)))*mphi * lyinv ,maxvalue)
        if (shear_real) then
          maxvalue=max(maxval(abs(grad_pot))*mphi * lyinv ,maxvalue)
        end if

        if (nlapar) then
          maxvalapar=max(maxval(abs(er(:,:,th)))*mrad*lxinv,maxvalapar)
          maxvalapar=max(maxval(abs(fr(:,:,th)))*mphi*lyinv,maxvalapar)
        end if

        if (nlbpar) then
          maxvalbpar=max(maxval(abs(gr(:,:,th)))*mrad*lxinv,maxvalbpar)
          maxvalbpar=max(maxval(abs(hr(:,:,th)))*mphi*lyinv,maxvalbpar)
        end if

        if (rhs_perf_measure) call perfoff(1)
      end if get_real_grad_phi

      ! Not until here is the parallel velocity looped over
      idxcopy2 = 1
      idx2 = 1

      !! Ideally need nested parallelism here
      !! Require runtime control of which threads go where.
      !! e.g. max possible into top level loop.  Then rest into this loop
      !! And two can be used for FFTs above.
      !! Idle threads will need to be utilised is section above
      !! The indexing arrays lindx3 and lincopy3 will need further dimensions
      !! indices.
      loop_nvpar: do kt = 1, nvpar

        if (rhs_perf_measure) call perfon('k-space multiplication 2',1)

        ! initialize a and b
        a(:,:,th) = (0.,0.)
        b(:,:,th) = (0.,0.)

        ! a = grad_x f_k
        ! b = grad_y f_k
        loop_nx2: do ix = 1, nx            ! Loop over radial modes
          cdum = ci1kxrh(ix)
          loop_nmod2: do imod = 1, nmod    ! Loop over poloidal modes
            !          a(imod,jind(ix),th) = ci1 * kxrh(ix)
            !                                    * fdis(indx(imod,ix,ipar,jv,kt,is))
            !          b(imod,jind(ix),th) = ci1 * krho(imod)
            !                                    * fdis(indx(imod,ix,ipar,jv,kt,is))
            !          a(imod,jind(ix),th) = ci1 * kxrh(ix) * fdis(lindx(idx))
            !          b(imod,jind(ix),th) = ci1 * krho(imod) * fdis(lindx(idx))
            !          idx=idx+1
            a(imod,jind(ix),th) =  cdum          * fdis(lindx3(idx2,ipar,jv,is))
            b(imod,jind(ix),th) =  ci1krho(imod) * fdis(lindx3(idx2,ipar,jv,is))
            idx2 = idx2 + 1
          end do loop_nmod2
        end do loop_nx2

        if (rhs_perf_measure) call perfoff(1)
        if (rhs_perf_measure) call perfon('fft 2',1)

        if(erase_any_transfer_from) then
          call erase_transfer_from_modes(a,th)
          call erase_transfer_from_modes(b,th)
        end if

        ! Inverse fourier transform k space to real space
        ! For gradients of distribution function
        if (non_linear) then
          call four2D_real(cr(:,:,th),a(:,:,th),FFT_INVERSE)
        else ! perp shear requires only the zeta component of the distribution
          cr(:,:,th) = 0.
        end if
        call four2D_real(dr(:,:,th),b(:,:,th),FFT_INVERSE)

        ! Now everything is in real space
        ! cr = grad_x f = radial gradient of the distribution in real space
        ! dr = grad_y f = zeta gradient of the distribution in real space

        if (rhs_perf_measure) call perfoff(1)
        if (rhs_perf_measure) call perfon('real-space multiplication',1)

        ! v_E grad f = (b x grad<phi> . grad f) /B
        ! As written in term III in the manual
        ! The factor of rhorat^2 comes from the normalisation of the FFT
        ! We neglect the d/ds variations
        ! Using also that the diagonal elements of efun tensor are 0:
        ! v_E grad f = (d<phi>/d(zeta))(df/d(psi))*
        !              E^zeta-phi+(d<phi>/d(psi))(df/d(zeta))E^phi-zeta

        ! The zeta psi component = efun(ipar,2,1)
        ! The psi zeta component = -efun(ipar,2,1) = efun(ipar,1,2)
        ! This routine is zero the flux tube version. Therefore use the 
        ! efun of the first radial grid point 
        dum = - efun(1,ipar,2,1)

        if (non_linear) then
          if (nlapar .and. (.not. nlbpar)) then
            do j = 1, mrad
              do i = 1, mphi
                ! cr = V_E . grad f = (b x grad <chi>).grad f
                ! The minus sign is due to the antisymmetry of efun
                ! efun is antisymmetric for all geometries, not just circular
                cr(i,j,th) = dum*((ar(i,j,th) &
                   &  - vpgr(ipar,jv,kt,is)*er(i,j,th))*cr(i,j,th) &
                   &  - ( br(i,j,th) &
                   &  - vpgr(ipar,jv,kt,is)*fr(i,j,th))*dr(i,j,th) )
              end do
            end do

          else if (nlbpar .and. (.not. nlapar)) then
            do j = 1, mrad
              do i = 1, mphi
                cr(i,j,th) = dum*((ar(i,j,th)+gr(i,j,th))*cr(i,j,th) &
                   & - (br(i,j,th)+hr(i,j,th))*dr(i,j,th))
              end do
            end do

          else if (nlapar .and. nlbpar) then
            do j = 1, mrad
              do i = 1, mphi
                cr(i,j,th) = dum*((ar(i,j,th)-vpgr(ipar,jv,kt,is)*er(i,j,th) &
                   & + gr(i,j,th))*cr(i,j,th) &
                   & - (br(i,j,th) - vpgr(ipar,jv,kt,is)*fr(i,j,th) &
                   & + hr(i,j,th))*dr(i,j,th))
              end do
            end do

          else
            do j = 1, mrad
              do i = 1, mphi
                ! cr = V_E . grad f = (b x grad <phi>).grad f
                ! The minus sign is due to the antisymmetry of efun
                ! efun is antisymmetric for all geometries, not just circular
                  cr(i,j,th) = dum*(ar(i,j,th)*cr(i,j,th)-br(i,j,th)*dr(i,j,th))
                ! ar = grad_zeta <phi> =
                ! poloidal gradient of the gyroaverage of potential in real space
                ! cr = grad_x f = radial gradient of the distribution in real space
                
                ! br = grad_psi <phi> =
                ! radial gradient of the gyroaverage of potential in real space
                ! dr = grad_y f = zeta gradient of the distribution in real space
              end do
            end do
          end if
        end if

        if (shear_real) then
          do j = 1, mrad
            do i = 1, mphi
              cr(i,j,th) = cr(i,j,th) + dum*grad_pot(j,ipar)*dr(i,j,th)
            end do
          end do
        end if

        if (rhs_perf_measure) call perfoff(1)
        if (rhs_perf_measure) call perfon('fft 3',1)

        ! Forward Fourier transform from real back to k space.
        ! Note that the positions of input and output arguments are reversed.
        call four2D_real(cr(:,:,th),a(:,:,th),FFT_FORWARD)

        ! Normalise - for more effciency do this operation after the copy so
        ! less elements are being normalised.. Also make division into a mult.
        a(:,:,th) = mphimrad1*a(:,:,th)
        ! Now the array a contains the sum of all nonlinear terms in the system,
        ! evaluated for this timestep and transformed to k space.

        if(erase_any_transfer_to) then
          call erase_transfer_to_modes(a,th)
        end if

        if (rhs_perf_measure) call perfoff(1)

        ! Copy the contribution into RHS
        ! The indexing arrays are for faster performance but make it 
        ! harder to see what is happening. 
        ! Practically the loops below do the following (identical but slower).
        ! The fact that here only the first (1:nmod,1:nx) modes
        ! are copied to the right-hand-side does the dealiasing.
        !do i = 1, nmod; do j = 1, nx
        !  if (jinv(j).ne.0) then
        !    if (.not.(jinv(j)==ixzero.and.i==1)) then
        !      rhs(indx(i,jinv(j),ipar,jv,kt,is)) =     & 
        !          &  rhs(indx(i,jinv(j),ipar,jv,kt,is)) + dtim*a(i,j,th) 
        !    end if
        !  endif
        !end do; end do

        if(present(rhs)) then
          ! Copy the result of the calculation into the rhs array:

          if (rhs_perf_measure) call perfon('copy',1)

          ! First start the i index from 2 to avoid adding in the (0,0) mode
          j = 1
          idxcopy2 = idxcopy2 + 1 ! start at 1 to avoid the first element
          do i = 2, nmod
            rhs(lincopy3(idxcopy2,ipar,jv,is)) =                                 &
               &    rhs(lincopy3(idxcopy2,ipar,jv,is)) + dtim_cmplx*a(i,j,th)
            idxcopy2 = idxcopy2 + 1
          end do

          do j = 2, nx-ixzero+1
            do i = 1, nmod
              rhs(lincopy3(idxcopy2,ipar,jv,is)) =                               &
                 &    rhs(lincopy3(idxcopy2,ipar,jv,is)) + dtim_cmplx*a(i,j,th)
              idxcopy2 = idxcopy2 + 1
            end do
          end do

          do j = mrad+2-ixzero, mrad
            do i = 1, nmod
              rhs(lincopy3(idxcopy2,ipar,jv,is)) =                               &
                 &    rhs(lincopy3(idxcopy2,ipar,jv,is)) + dtim_cmplx*a(i,j,th)
              idxcopy2 = idxcopy2 + 1
            end do
          end do

          if (rhs_perf_measure) call perfoff(1)
        end if

    end do loop_nvpar

   ! end do loop_ns
   ! end do loop_nmu
   ! end do loop_nsp
 end do loop3

  
!$OMP end parallel

if(present(rhs)) then
  if (rhs_perf_measure) call perfon('Estimate timestep limit',1)

  ! \attention maxvalue can be zero at this point (and is so in some of the
  !            testcases).
  if(abs(maxvalue) >= r_tiny) then
    dtim_est = 2./maxvalue
  else
    dtim_est = r_huge
  end if
  ! Save minimum timestep so far for local processor
  dtim_est_save = min(dtim_est,dtim_est_save)

  if (number_of_processors > 1 .and. nl_dtim_est) then
    call mpiallreduce_min(dtim_est,dtim_est_dum,1)
    dtim_est = dtim_est_dum
  end if

  if (nlapar) then

    ! rather crude way to measure the time step limit for electromagnetic
    ! cases
    maxvalapar = maxvalapar*vpmax

    if (maxvalapar > tiny(1.E0)) dtim_est_apar = 2./maxvalapar

    ! Save minimum timestep so far for local processor
    dtim_est_save = min(dtim_est_apar,dtim_est_save)

    if (number_of_processors > 1 .and. nl_dtim_est) then
      call mpiallreduce_min(dtim_est_apar,dtim_est_dum,1)
      dtim_est = min(dtim_est,dtim_est_dum)
    end if

  end if
  if (rhs_perf_measure) call perfoff(1)
end if

if (rhs_perf_measure) call perfoff(1)

end subroutine add_non_linear_terms_spectral



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Input: a field in k-space representation.
!> This subroutine sets the modes to zero which specified in
!> mode::no_transfer_from_modes
!-------------------------------------------------------------
pure subroutine erase_transfer_from_modes(a, th)
  use mode, only : no_transfer_from_modes
  complex, intent(inout) :: a(:,:,:)
  integer, intent(in) :: th
  integer :: i
  do i = 1, size(no_transfer_from_modes,2)
    if(no_transfer_from_modes(1,i) == 0) then
      a(:, jind(no_transfer_from_modes(2,i)), th) = 0
    elseif(no_transfer_from_modes(2,i) == 0) then
      a(no_transfer_from_modes(1,i), :, th) = 0
    else
      a(no_transfer_from_modes(1,i), jind(no_transfer_from_modes(2,i)), th) = 0
    end if
  end do
  
end subroutine erase_transfer_from_modes


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Input: a field in k-space representation.
!> This subroutine sets the modes to zero which specified in
!> mode::no_transfer_to_modes
!-------------------------------------------------------------
pure subroutine erase_transfer_to_modes(a, th)
  use mode, only : no_transfer_to_modes
  complex, intent(inout) :: a(:,:,:)
  integer, intent(in) :: th
  integer :: i
  do i = 1, size(no_transfer_to_modes,2)
    if(no_transfer_to_modes(1,i) == 0) then
      a(:, jind(no_transfer_to_modes(2,i)), th) = 0
    elseif(no_transfer_to_modes(2,i) == 0) then
      a(no_transfer_to_modes(1,i), :, th) = 0
    else
      a(no_transfer_to_modes(1,i), jind(no_transfer_to_modes(2,i)), th) = 0
    end if
  end do
  
end subroutine erase_transfer_to_modes

!------------------------------------------------------------------------------
!> This routine adds the nonlinear terms for the non_spectral case. 
!> This routine uses the Arakawa scheme (fourth order). Potential and 
!> distribution function are first Fourier transformed to real space where
!> the Arakawa scheme is applied and then a back transformation is 
!> performed. 
!------------------------------------------------------------------------------
subroutine add_non_linear_terms_nonspec(fdis,rhs) 

  use control,        only : non_linear, nlapar, nlbpar, dtim, dtim_est
  use control,        only : nl_dtim_est, fac_dtim_nl, dtim_est_save
  use control,        only : dtim
  use dist,           only : nsolc, msolc
  use matdat,         only : get_f_from_g
  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp, lsendrecv_x, vpmax
  use grid,           only : lx, n_x_grid 
  use ompinterface,   only : ompget_thread_num
  use mpiinterface,   only : number_of_processors, mpiallreduce_min
  use general,        only : gkw_abort
  use fft,            only : four1D_real, FFT_INVERSE, FFT_FORWARD
  use geom,           only : efun, dxgr
  use mode,           only : lyinv, lxinv, krho
  use mode,           only : erase_any_transfer_to, erase_any_transfer_from
  use components,     only : vthrat, signz
  use velocitygrid,   only : vpgr
  use global,         only : r_tiny, r_huge
  use constants,      only : c1, pi

  complex, intent(in)    :: fdis(:)
  complex, optional, intent(inout) :: rhs(nsolc)

  integer :: imod, ix, i, j, k, is, th, im, ir, idx, idx2, iloop, n_fields
  integer :: ifields 
  complex :: dtim_cmplx
  real    :: dum, dom, maxvalue, dtim_est_dum, maxvalapar
  !The maximum aparallel per tordoidal mode for the time step estimator.
  real    :: maxapar_pm(nmod)

  ! return if called for shear_real (this is added in linear terms)
  if (.not. non_linear) return

  if (.not. nl_initialised) then
    call gkw_abort('add_non_linear terms: cannot call before nonlinear_init')
  end if

  ! checks on the size of passed fdis (could be removed)
  if (lsendrecv_x) then
    if (size(fdis) /= msolc) call gkw_abort('nl term wrong size fdis') 
  else
    if (size(fdis) /= nsolc) call gkw_abort('nl term wrong fdis size') 
  end if

  ! stuff that does not work at present (nlbpar) 
  if (nlbpar) then 
    call gkw_abort('For the nonspectral method nlbpar is not implemented')
  end if 

  ! set the value for the estimate of the timestep to zero (The restriction
  ! due to the electromagnetic scheme are measured separately with maxvalapar)
  maxapar_pm = 0.
  maxvalue   = 0.
  maxvalapar = 0.
 
  ! use complex dtim
  dtim_cmplx = c1*dtim

  n_fields = 1
  if (nlapar) n_fields = 2

  ! outer loop over i, j, is 
  !do is = 1, nsp; do j = 1, nmu; do i = 1, ns  
 
  ! Manually collapse the three nested loops above
  !!$OMP do collapse(3) Not supported until OpenMP 3.0
 
  !$OMP parallel do private(iloop,i,is,j,ix,imod,k,th,idx2,idx,im, &
  !$OMP ir)
  loop3: do iloop = 1, ns*nsp*nmu

    is = i3loop(1,iloop)
    j = i3loop(2,iloop)
    i = i3loop(3,iloop)

    ! Always 1 for non openmp runs
    th = ompget_thread_num() + 1

    idx=1

      ! fill the fields arrays
      do ifields = n_fields,1,-1 ! the linphi array returns iapar_ga first

        ! fill the array a with the potential 
        a(:,:,th) = (0.E0,0.E0) 

        do ix = 1-xgp, nx+xgp
          do imod = 1, nmod
            a(imod,ix,th) = fdis(linphi(idx,i,j,is))
            idx = idx + 1
          end do
        end do

        if(ifields.eq.2)then
          do imod = 1,nmod
            maxapar_pm(imod) = max(krho(imod)*maxval(abs(a(imod,:,th))),& 
              & maxapar_pm(imod))
          end do
          maxvalapar = max(maxval(maxapar_pm(:)),maxvalapar)
        endif

        ! the boundary condition
        do imod = 1, nmod
          do ix = 1, x_boundary_points
            a(imod,  1-ix,th) = mbd(1-ix,imod,i)*a(imod,ibd(1-ix,imod,i),th) 
            a(imod,nx+ix,th) = mbd(ix,imod,i)*a(imod,ibd( ix,imod,i),th)
          end do
        end do

        if(erase_any_transfer_from) then
          call erase_transfer_from_modes(a,th)        
        end if


        ! do the inverse Fourier transform 
        ! FJC: Faster here to use fourcol -> fourcol_ra_ca ?
        ! which does all the rows of the array in 1D together
        do ix = 1-x_boundary_points, nx+x_boundary_points  
          call four1D_real(ar(1:mphi,ix,th),a(:,ix,th),FFT_INVERSE)
        end do
         
        ! Minimal time saving (~ 2%) and not working
        !call fourcol(ar(1:mphi,1-x_boundary_points:nx+x_boundary_points,th),&
        !     &a(:,1-x_boundary_points:nx+x_boundary_points,th),FFT_INVERSE)

        ! the periodic boundary conditions in the zeta direction 
        do ix = 1-x_boundary_points, nx+x_boundary_points 
          ar(    -1,ix,th) = ar(mphi-1,ix,th) 
          ar(     0,ix,th) = ar(  mphi,ix,th)
          ar(mphi+1,ix,th) = ar(     1,ix,th)
          ar(mphi+2,ix,th) = ar(     2,ix,th) 
        end do

        ! store iapar_ga in er (this happens first time
        if (ifields == 2) er(:,:,th) =  ar(:,:,th)
        ! For electromagnetic runs the potential needs to be stored separately
        if (nlapar) then
          if (ifields == 1) fr(:,:,th) = ar(:,:,th) 
        endif
        

      end do ! ifields loop

      ! loop over the parallel velocity
      idx=1; idx2=1
      vpar_loop : do k = 1, nvpar

        ! copy out the distribution function 
        b(:,:,th) = (0.E0,0.E0)
        do ix = 1-xgp, nx+xgp; do imod = 1, nmod
          !b(imod,ix,th) = fdis(indx(ifdis,imod,ix,i,j,k,is))
          b(imod,ix,th) = fdis(lindx3(idx,i,j,is))
          idx = idx + 1
        end do; end do; 

        ! the boundary conditions 
        do imod = 1, nmod
          do ix = 1, x_boundary_points
            b(imod,  1-ix,th) = mbd(1-ix,imod,i)*b(imod,ibd(1-ix,imod,i),th)
            b(imod,nx+ix,th) = mbd(ix,imod,i)*b(imod,ibd(ix,imod,i),th)
          end do
        end do

        if(erase_any_transfer_from) then
          call erase_transfer_from_modes(b,th)
        end if

        ! do the Fourier transform to position-space
        do ix = 1-x_boundary_points, nx+x_boundary_points
          call four1D_real(br(1:mphi,ix,th),b(:,ix,th),FFT_INVERSE)
        end do

        !call fourcol(br(1:mphi,1-x_boundary_points:nx+x_boundary_points,th),&
        !     &b(:,1-x_boundary_points:nx+x_boundary_points,th),FFT_INVERSE)

        ! the periodic boundary conditions in the zeta direction 
        do ix = 1-x_boundary_points, nx+x_boundary_points 
          br(    -1,ix,th) = br(mphi-1,ix,th) 
          br(     0,ix,th) = br(  mphi,ix,th)
          br(mphi+1,ix,th) = br(     1,ix,th)
          br(mphi+2,ix,th) = br(     2,ix,th) 
        end do

        ! put chi combined field into ar for electromagnetic runs 
        if (nlapar)  ar(:,:,th) = fr(:,:,th) - 2. * vpgr(i,j,k,is) * vthrat(is) * er(:,:,th)

        if (.true.) then 
          ! the first bracket
          if (.not. nl_dtim_est) then  
            do ir = 1, nx ; do im = 1, mphi; 
              cr(im,ir,th)= &
                 (ar(im-2,ir,th)-8.0E0*ar(im-1,ir,th)+8.0E0*ar(im+1,ir,th)-ar(im+2,ir,th))* &
                 (br(im,ir-2,th)-8.0E0*br(im,ir-1,th)+8.0E0*br(im,ir+1,th)-br(im,ir+2,th))- &
                 (br(im-2,ir,th)-8.0E0*br(im-1,ir,th)+8.0E0*br(im+1,ir,th)-br(im+2,ir,th))* &
                 (ar(im,ir-2,th)-8.0E0*ar(im,ir-1,th)+8.0E0*ar(im,ir+1,th)-ar(im,ir+2,th))
            end do; end do   
          else
            do ir = 1, nx ; do im = 1, mphi; 
              dum =  ar(im-2,ir,th)-8.0E0*ar(im-1,ir,th)+8.0E0*ar(im+1,ir,th)-ar(im+2,ir,th)
              dom =  ar(im,ir-2,th)-8.0E0*ar(im,ir-1,th)+8.0E0*ar(im,ir+1,th)-ar(im,ir+2,th)
              cr(im,ir,th) = dum * &
                 (br(im,ir-2,th)-8.0E0*br(im,ir-1,th)+8.0E0*br(im,ir+1,th)-br(im,ir+2,th))- &
                 dom * &
                 (br(im-2,ir,th)-8.0E0*br(im-1,ir,th)+8.0E0*br(im+1,ir,th)-br(im+2,ir,th))  
              maxvalue=MAX(abs(dom)*mrad*lxinv/12.E0,maxvalue)
              maxvalue=MAX(abs(dum)*mphi*lyinv/12.E0,maxvalue)
            end do; end do   
          endif

          ! second bracket 
          do ir = 1, nx;  do im = -1, mphi+2 
            dr(im,ir,th) = ar(im,ir,th) * &
               (br(im,ir-2,th)-8.0E0*br(im,ir-1,th)+8.0E0*br(im,ir+1,th)-br(im,ir+2,th))
          end do; end do; 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) + &  
               (dr(im-2,ir,th)-8.0E0*dr(im-1,ir,th)+8.0E0*dr(im+1,ir,th)-dr(im+2,ir,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi 
            dr(im,ir,th) = ar(im,ir,th) * & 
               (br(im-2,ir,th)-8.0E0*br(im-1,ir,th)+8.0E0*br(im+1,ir,th)-br(im+2,ir,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) - &  
               (dr(im,ir-2,th)-8.0E0*dr(im,ir-1,th)+8.0E0*dr(im,ir+1,th)-dr(im,ir+2,th))
          end do; end do 

          ! third bracket 
          do ir = 1, nx; do im = -1, mphi+2
            dr(im,ir,th) = br(im,ir,th) * &
               (ar(im,ir-2,th)-8.0E0*ar(im,ir-1,th)+8.0E0*ar(im,ir+1,th)-ar(im,ir+2,th))
          end do; end do; 
          do ir = 1, nx; do im = 1, mphi 
            cr(im,ir,th) = cr(im,ir,th) - &  
               (dr(im-2,ir,th)-8.0E0*dr(im-1,ir,th)+8.0E0*dr(im+1,ir,th)-dr(im+2,ir,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi
            dr(im,ir,th) = br(im,ir,th) * & 
               (ar(im-2,ir,th)-8.0E0*ar(im-1,ir,th)+8.0E0*ar(im+1,ir,th)-ar(im+2,ir,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) + &  
               (dr(im,ir-2,th)-8.0E0*dr(im,ir-1,th)+8.0E0*dr(im,ir+1,th)-dr(im,ir+2,th))
          end do; end do 

        else ! second order 

          ! the first bracket 
          do ir = 1, nx ; do im = 1, mphi; 
            cr(im,ir,th)= &
               (-6.0E0*ar(im-1,ir,th)+6.0E0*ar(im+1,ir,th))* &
               (-6.0E0*br(im,ir-1,th)+6.0E0*br(im,ir+1,th))- &
               (-6.0E0*br(im-1,ir,th)+6.0E0*br(im+1,ir,th))* &
               (-6.0E0*ar(im,ir-1,th)+6.0E0*ar(im,ir+1,th))
          end do; end do   

          ! second bracket 
          do ir = 1, nx;  do im = -1, mphi+2 
            dr(im,ir,th) = ar(im,ir,th) * &
               (-6.0E0*br(im,ir-1,th)+6.0E0*br(im,ir+1,th))
          end do; end do; 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) + &  
               (-6.0E0*dr(im-1,ir,th)+6.0E0*dr(im+1,ir,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi 
            dr(im,ir,th) = ar(im,ir,th) * & 
               (-6.0E0*br(im-1,ir,th)+6.0E0*br(im+1,ir,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) - &  
               (-6.0E0*dr(im,ir-1,th)+6.0E0*dr(im,ir+1,th))
          end do; end do 

          ! third bracket 
          do ir = 1, nx; do im = -1, mphi+2
            dr(im,ir,th) = br(im,ir,th) * &
               (-6.0E0*ar(im,ir-1,th)+6.0E0*ar(im,ir+1,th))
          end do; end do; 
          do ir = 1, nx; do im = 1, mphi 
            cr(im,ir,th) = cr(im,ir,th) - &  
               (-6.0E0*dr(im-1,ir,th)+6.0E0*dr(im+1,ir,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi
            dr(im,ir,th) = br(im,ir,th) * & 
               (-6.0E0*ar(im-1,ir,th)+6.0E0*ar(im+1,ir,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi
            cr(im,ir,th) = cr(im,ir,th) + &  
               (-6.0E0*dr(im,ir-1,th)+6.0E0*dr(im,ir+1,th))
          end do; end do 

        endif

        ! do the back transform to k-space
        do ix = 1, nx 
          call four1D_real(cr(1:mphi,ix,th),a(:,ix,th),FFT_FORWARD)
        end do

        if(erase_any_transfer_to) then
          call erase_transfer_to_modes(a,th)
        end if

        ! cannot use this as the array is not contiguous in memory
        !call fourcol(cr(1:mphi,1:nx,th),a(:,1:nx,th),FFT_FORWARD)

        ! get the size of 0,0 mode (should be tiny)
        ! reducing this over x could be expensive unless done outside loop
        ! not yet demonstrated that this is needed
        !dumc=sum(a(1,1:nx,th))/real(nx)
        !dumc=(0.,0.)
        if(present(rhs)) then
          ! add to rhs (may need to remove the 0,0 average)
          do ix = 1, nx; do imod = 1, nmod 
            rhs(lincopy3(idx2,i,j,is)) = rhs(lincopy3(idx2,i,j,is)) - &
               dtim_cmplx * efun(ix,i,2,1) * lyinv * a(imod,ix,th)  / (432.0E0 *  &
               dxgr )
            idx2 = idx2 + 1
          end do; end do;
        end if

      end do vpar_loop

    end do loop3
    !$OMP end parallel do

    if(present(rhs)) then
  ! \attention maxvalue can be zero at this point (and will be in some of the
  !            testcases).
  if(abs(maxvalue) >= r_tiny) then
    dtim_est = 2./maxvalue
  else
    dtim_est = r_huge
  end if

  !This is the nonlinear timestep estimator for electromagnetic
  !global runs.
  !Rhostar doesnt need to be here explicitly as it is contained
  !within lx 

  if(nlapar.and.nl_dtim_est)then
    do is = 1,nsp
      !We should only consider electrons here
      if(signz(is).lt.0)then
        dum  = 2.E0*pi*vpmax*vthrat(is)*n_x_grid/lx 
        !Maybe there should be a Bn here too.
        dtim_est = min(dtim_est,1.0/(maxvalapar*dum))
      end if
    end do
  end if

  ! Save minimum timestep so far for local processor
  dtim_est_save = min(dtim_est,dtim_est_save)

  if (number_of_processors > 1 .and. nl_dtim_est) then
    call mpiallreduce_min(dtim_est_save,dtim_est_dum,1)
    dtim_est = dtim_est_dum
  end if

  dtim_est=dtim_est*fac_dtim_nl
end if

end subroutine add_non_linear_terms_nonspec


!------------------------------------------------------------------------------
!> This routine adds the nonlinear terms for the non_spectral case. 
!> This routine uses the Arakawa scheme (fourth order) potential and 
!> distribution function are first Fourier transformed to real space 
!> the Arakawa scheme is applied and then a back transformation is 
!> performed. 
!> This is the second of the nonspectral routines. It includes the Sung 
!> effect and the velocity nonlinearity. It has been optimised 
!------------------------------------------------------------------------------
subroutine add_non_linear_terms_nonspec_2(fdis,rhs) 

  use components,     only : vthrat, rhostar
  use control,        only : non_linear, nlapar, nlbpar, dtim, dtim_est
  use control,        only : nl_dtim_est, lpar_vel_nl
  use dist,           only : nsolc, msolc
  use fft,            only : four1D_real, FFT_INVERSE, FFT_FORWARD
  use general,        only : gkw_abort, gkw_warn
  use geom,           only : efun, dxgr, sgr_dist
  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp
  use grid,           only : lsendrecv_x 
  use mode,           only : lyinv
  use mode,           only : erase_any_transfer_to, erase_any_transfer_from
  use mode,           only : no_transfer_to_modes, no_transfer_from_modes
  use ompinterface,   only : ompget_thread_num
  use velocitygrid,   only : vpgr
  use rho_par_switch, only : lnonlinear_rhostar 

  complex, intent(in)    :: fdis(:)
  complex, intent(inout) :: rhs(nsolc)

  integer :: imod, ix, i, j, k, l, is, th, im, ir, idx, idx2, iloop, n_fields
  integer :: ifields, inp, inpf, ipe

  ! return if called for shear_real (this is added in linear terms)
  if (.not. non_linear) return

  if (.not. nl_initialised) then
    call gkw_abort('add_non_linear terms: cannot call before nonlinear_init')
  end if

  ! checks on the size of passed fdis (could be removed)
  if (lsendrecv_x) then
    if (size(fdis) /= msolc) call gkw_abort('nl term wrong size fdis') 
  else
    if (size(fdis) /= nsolc) call gkw_abort('nl term wrong fdis size') 
  end if

  ! stuff that does not work at pressent (nlbpar) 
  if (nlbpar) then 
    call gkw_abort('For the nonspectral method nlbpar is not implemented')
  end if 
  
  if (nl_dtim_est) then 
    call gkw_warn('At present no estimate of NL timestep in nonspectral')
    nl_dtim_est = .false.
    dtim_est = dtim + 1.0
  end if 

  n_fields = 1
  if (nlapar) n_fields = 2
 
  ! Manually collapse the three nested loops above
  !!$OMP do collapse(3) Not supported until OpenMP 3.0
 
  !!!$OMP parallel do default(shared) private(iloop,i,is,j,ix,imod,k,th,idx2,idx,im,ir,dumc)
  
  ! set the counter for the parallel velocity nonlinearity 
  inpf = 0 
  ipe  = 0 

  loop3: do iloop = 1, nsp*nmu

    is = i3loop(1,iloop)
    j  = i3loop(2,iloop)

    ! Always 1 for non openmp runs
    th = ompget_thread_num() + 1
    idx=1

    ! fill the fields arrays
    do ifields = n_fields,1,-1 ! the linphi array returns iapar_ga first

      ! fill the array a with the field  
      aa(:,:,:,th) = (0.E0,0.E0) 

      ! outer loop over the parallel grid points 
      do i = 1, ns 

        do ix = 1-xgp, nx+xgp; do imod = 1, nmod   
          aa(imod,ix,i,th) = fdis(linphi(idx,1,j,is))
          idx = idx + 1
        end do; end do      

        ! the boundary condition
        do imod = 1, nmod 
          aa(imod,  -1,i,th) = mbd(-1,imod,i)*aa(imod,ibd(-1,imod,i),i,th) 
          aa(imod,   0,i,th) = mbd( 0,imod,i)*aa(imod,ibd( 0,imod,i),i,th) 
          aa(imod,nx+1,i,th) = mbd( 1,imod,i)*aa(imod,ibd( 1,imod,i),i,th) 
          aa(imod,nx+2,i,th) = mbd( 2,imod,i)*aa(imod,ibd( 2,imod,i),i,th) 
        end do 

        if(erase_any_transfer_from) then
          do l = 1, size(no_transfer_from_modes)
            aa(no_transfer_from_modes(1,i),:, :, th) = 0
          end do
        end if
        
        ! do the inverse Fourier transform 
        do ix = 1-x_boundary_points, nx+x_boundary_points   
          call four1D_real(aar(1:mphi,ix,i,th),aa(:,ix,i,th),FFT_INVERSE)
        end do 
      
        ! the periodic boundary conditions in the zeta direction 
        do ix = 1-x_boundary_points, nx+x_boundary_points  
          aar(    -1,ix,i,th) = aar(mphi-1,ix,i,th) 
          aar(     0,ix,i,th) = aar(  mphi,ix,i,th)
          aar(mphi+1,ix,i,th) = aar(     1,ix,i,th)
          aar(mphi+2,ix,i,th) = aar(     2,ix,i,th) 
        end do 

      end do ! loop over parallel grid points 

      ! store iapar_ga in er (this happens first time
      if (ifields == 2) eer(:,:,:,th) =  aar(:,:,:,th)
      ! For electromagnetic runs the potential needs to be stored separately
      if (nlapar) then 
        if (ifields == 1) ffr(:,:,:,th) = aar(:,:,:,th) 
      endif 

    end do ! ifields loop

    ! loop over the parallel velocity
    idx=1; idx2=1
    vpar_loop : do k = 1, nvpar
     
      ! copy out the distribution function 
      bb(:,:,:,th) = (0.E0,0.E0)

      ! outer loop over the parallel grid points 
      do i = 1, ns 

        do ix = 1-xgp, nx+xgp; do imod = 1, nmod 
          bb(imod,ix,i,th) = fdis(lindx3(idx,1,j,is))
          idx = idx + 1
        end do; end do  

        ! the boundary conditions 
        do imod = 1, nmod 
          bb(imod,  -1,i,th) = mbd(-1,imod,i)*bb(imod,ibd(-1,imod,i),i,th) 
          bb(imod,   0,i,th) = mbd( 0,imod,i)*bb(imod,ibd( 0,imod,i),i,th) 
          bb(imod,nx+1,i,th) = mbd( 1,imod,i)*bb(imod,ibd( 1,imod,i),i,th) 
          bb(imod,nx+2,i,th) = mbd( 2,imod,i)*bb(imod,ibd( 2,imod,i),i,th) 
        end do 

        if(erase_any_transfer_from) then
          do l = 1, size(no_transfer_from_modes)
            bb(no_transfer_from_modes(1,i),:, :, th) = 0
          end do
        end if
        
        ! do the inverse Fourier transform 
        do ix = 1-x_boundary_points, nx+x_boundary_points
          call four1D_real(bbr(1:mphi,ix,i,th),bb(:,ix,i,th),FFT_INVERSE)
        end do 
        
        ! the periodic boundary conditions in the zeta direction 
        do ix = 1-x_boundary_points, nx+x_boundary_points 
          bbr(    -1,ix,i,th) = bbr(mphi-1,ix,i,th) 
          bbr(     0,ix,i,th) = bbr(  mphi,ix,i,th)
          bbr(mphi+1,ix,i,th) = bbr(     1,ix,i,th)
          bbr(mphi+2,ix,i,th) = bbr(     2,ix,i,th) 
        end do 

        ! put chi combined field into ar for electromagnetic runs 
        if (nlapar)  then 
            aar(:,:,i,th) = ffr(:,:,i,th) - 2. * vpgr(i,j,k,is) *  & 
                          & vthrat(is) * eer(:,:,i,th)
        end if 
  
        tt: if (.true.) then 
          ! the first bracket 
          do ir = 1, nx ; do im = 1, mphi   
            ccr(im,ir,i,th)= &
              (aar(im-2,ir,i,th)-8.0E0*aar(im-1,ir,i,th)+8.0E0*aar(im+1,ir,i,th)-aar(im+2,ir,i,th))* &
              (bbr(im,ir-2,i,th)-8.0E0*bbr(im,ir-1,i,th)+8.0E0*bbr(im,ir+1,i,th)-bbr(im,ir+2,i,th))- &
              (bbr(im-2,ir,i,th)-8.0E0*bbr(im-1,ir,i,th)+8.0E0*bbr(im+1,ir,i,th)-bbr(im+2,ir,i,th))* &
              (aar(im,ir-2,i,th)-8.0E0*aar(im,ir-1,i,th)+8.0E0*aar(im,ir+1,i,th)-aar(im,ir+2,i,th))
          end do; end do    

          ! second bracket 
          do ir = 1, nx;  do im = -1, mphi+2  
            ddr(im,ir,i,th) = aar(im,ir,i,th) * &
              (bbr(im,ir-2,i,th)-8.0E0*bbr(im,ir-1,i,th)+8.0E0*bbr(im,ir+1,i,th)-bbr(im,ir+2,i,th))
          end do; end do; 
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) + &  
              (ddr(im-2,ir,i,th)-8.0E0*ddr(im-1,ir,i,th)+8.0E0*ddr(im+1,ir,i,th)-ddr(im+2,ir,i,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi  
            ddr(im,ir,i,th) = aar(im,ir,i,th) * & 
              (bbr(im-2,ir,i,th)-8.0E0*bbr(im-1,ir,i,th)+8.0E0*bbr(im+1,ir,i,th)-bbr(im+2,ir,i,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) - &  
              (ddr(im,ir-2,i,th)-8.0E0*ddr(im,ir-1,i,th)+8.0E0*ddr(im,ir+1,i,th)-ddr(im,ir+2,i,th))
          end do; end do 

          ! third bracket 
          do ir = 1, nx; do im = -1, mphi+2 
            ddr(im,ir,i,th) = bbr(im,ir,i,th) * &
              (aar(im,ir-2,i,th)-8.0E0*aar(im,ir-1,i,th)+8.0E0*aar(im,ir+1,i,th)-aar(im,ir+2,i,th))
          end do; end do 
          do ir = 1, nx; do im = 1, mphi 
            ccr(im,ir,i,th) = ccr(im,ir,i,th) - &  
              (ddr(im-2,ir,i,th)-8.0E0*ddr(im-1,ir,i,th)+8.0E0*ddr(im+1,ir,i,th)-ddr(im+2,ir,i,th))
          end do; end do 
          do ir = -1, nx+2; do im = 1, mphi 
            ddr(im,ir,i,th) = bbr(im,ir,i,th) * & 
              (aar(im-2,ir,i,th)-8.0E0*aar(im-1,ir,i,th)+8.0E0*aar(im+1,ir,i,th)-aar(im+2,ir,i,th))
          end do; end do  
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) + &  
              (ddr(im,ir-2,i,th)-8.0E0*ddr(im,ir-1,i,th)+8.0E0*ddr(im,ir+1,i,th)-ddr(im,ir+2,i,th))
          end do; end do  

        else ! second order 

          ! the first bracket 
          do ir = 1, nx ; do im = 1, mphi  
            ccr(im,ir,i,th)= &
              (-6.0E0*aar(im-1,ir,i,th)+6.0E0*aar(im+1,ir,i,th))* &
              (-6.0E0*bbr(im,ir-1,i,th)+6.0E0*bbr(im,ir+1,i,th))- &
              (-6.0E0*bbr(im-1,ir,i,th)+6.0E0*bbr(im+1,ir,i,th))* &
              (-6.0E0*aar(im,ir-1,i,th)+6.0E0*aar(im,ir+1,i,th))
          end do; end do 

          ! second bracket 
          do ir = 1, nx;  do im = -1, mphi+2 
            ddr(im,ir,i,th) = aar(im,ir,i,th) * &
              (-6.0E0*bbr(im,ir-1,i,th)+6.0E0*bbr(im,ir+1,i,th))
          end do; end do  
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) + &  
              (-6.0E0*ddr(im-1,ir,i,th)+6.0E0*ddr(im+1,ir,i,th))
          end do; end do  
          do ir = -1, nx+2; do im = 1, mphi
            ddr(im,ir,i,th) = aar(im,ir,i,th) * & 
              (-6.0E0*bbr(im-1,ir,i,th)+6.0E0*bbr(im+1,ir,i,th))
          end do; end do  
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) - &  
              (-6.0E0*ddr(im,ir-1,i,th)+6.0E0*ddr(im,ir+1,i,th))
          end do; end do  

          ! third bracket 
          do ir = 1, nx; do im = -1, mphi+2 
            ddr(im,ir,i,th) = bbr(im,ir,i,th) * &
              (-6.0E0*aar(im,ir-1,i,th)+6.0E0*aar(im,ir+1,i,th))
          end do; end do  
          do ir = 1, nx; do im = 1, mphi 
            ccr(im,ir,i,th) = ccr(im,ir,i,th) - &  
              (-6.0E0*ddr(im-1,ir,i,th)+6.0E0*ddr(im+1,ir,i,th))
          end do; end do  
          do ir = -1, nx+2; do im = 1, mphi 
            ddr(im,ir,i,th) = bbr(im,ir,i,th) * & 
              (-6.0E0*aar(im-1,ir,i,th)+6.0E0*aar(im+1,ir,i,th))
          end do; end do  
          do ir = 1, nx; do im = 1, mphi
            ccr(im,ir,i,th) = ccr(im,ir,i,th) + &  
              (-6.0E0*ddr(im,ir-1,i,th)+6.0E0*ddr(im,ir+1,i,th))
          end do; end do 

        endif tt

        do ir = 1, nx; do im = 1, mphi
          ccr(im,ir,i,th) = ccr(im,ir,i,th) * dtim * efun(ir,i,2,1) * lyinv / (432.0E0 * dxgr )
        end do; end do 

      end do ! loop over the parallel grid points

      
      pvnl: if (lpar_vel_nl) then 
         
        ! loop over the parallel grid points 
        do i = 1, ns       

          ! fill the bb array with the velocity derivative of f  
          do ix = 1, nx; do imod = 1, nmod 
            inpf = inpf + 1           
            bb(imod,ix,i,th) = cpvnl(1,inpf)*fdis(ipvnl(1,inpf))+cpvnl(2,inpf)*fdis(ipvnl(2,inpf))+  &
                             & cpvnl(3,inpf)*fdis(ipvnl(3,inpf))+cpvnl(4,inpf)*fdis(ipvnl(4,inpf))
          end do; end do  
      
          ! do the inverse Fourier transform 
          do ix = 1, nx
            call four1D_real(bbr(1:mphi,ix,i,th),bb(:,ix,i,th),FFT_INVERSE)
          end do 
      
        end do 
   
        ! set the counter 
        inp = 0 
        
        ! calculate the inner product of dX/dt and nabla \phi
        ! First the parallel derivative 
        do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
            inp = inp + 1 
            ddr(imod,ix,i,th) =  fpvnl(is) *                              &
                              &  real(dpvnl(1,inp)*aar(imod,ix,ipv(1,i),th)   &
                              &  +dpvnl(2,inp)*aar(imod,ix,ipv(2,i),th)   &
                              &  +dpvnl(3,inp)*aar(imod,ix,ipv(3,i),th)   &
                              &  +dpvnl(4,inp)*aar(imod,ix,ipv(4,i),th) ) 
        end do; end do; end do 
        
        ! the drift terms 
        do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
          ipe = ipe + 1
          ! aar is backtrafo of fouriertrafo of a real valued function and
          ! therefore should also be real valued. The cast to real is thus
          ! assumed to be safe.
          ddr(imod,ix,i,th) =  ddr(imod,ix,i,th) +                         &
                            &  epvnl(1,ipe)*real(    aar(imod-2,ix,i,th)   &
                            &                - 8.0E0*aar(imod-1,ix,i,th)   & 
                            &                + 8.0E0*aar(imod+1,ix,i,th)   &
                            &                      - aar(imod+2,ix,i,th) ) &
                            & +epvnl(2,ipe)*real(    aar(imod,ix-2,i,th)   &
                            &                - 8.0E0*aar(imod,ix-1,i,th)   & 
                            &                + 8.0E0*aar(imod,ix+1,i,th)   & 
                            &                      - aar(imod,ix+2,i,th) )
        end do; end do; end do 

        ! now add the velocity nonlinearity to ccr 
        do i = 1, ns; do ix = 1, nx; do imod = 1, nmod 
          ccr(imod,ix,i,th) = ccr(imod,ix,i,th) + dtim * ddr(imod,ix,i,th) &
                            & *bbr(imod,ix,i,th) 
        end do; end do; end do         
        
      endif pvnl 
      
      nlr: if (lnonlinear_rhostar) then 

        ! the first bracket 
        do i = 1, ns; do ir = 1, nx ; do im = 1, mphi   
          ggr(im,ir,i,th)= &
            real(cps(i,-2)*aar(im,ir,ips(i,-2),th)+cps(i,-1)*aar(im,ir,ips(i,-1),th)+              &
            & cps(i,1)*aar(im,ir,ips(i,1),th)+cps(i,2)*aar(im,ir,ips(i,2),th)   ) *                &
            (bbr(im,ir-2,i,th)-8.0E0*bbr(im,ir-1,i,th)+8.0E0*bbr(im,ir+1,i,th)-bbr(im,ir+2,i,th))- &
            real(cps(i,-2)*bbr(im,ir,ips(i,-2),th)+cps(i,-1)*bbr(im,ir,ips(i,-1),th)+              &
            & cps(i,1)*bbr(im,ir,ips(i,1),th)+cps(i,2)*bbr(im,ir,ips(i,2),th)   ) *                &
            (aar(im,ir-2,i,th)-8.0E0*aar(im,ir-1,i,th)+8.0E0*aar(im,ir+1,i,th)-aar(im,ir+2,i,th))
        end do; end do; end do    

        ! second bracket 
        do i = 1, ns; do ir = 1, nx;  do im = -1, mphi+2  
          ddr(im,ir,i,th) = aar(im,ir,i,th) * &
            (bbr(im,ir-2,i,th)-8.0E0*bbr(im,ir-1,i,th)+8.0E0*bbr(im,ir+1,i,th)-bbr(im,ir+2,i,th))
        end do; end do; end do; 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) + &  
            real(cps(i,-2)*ddr(im,ir,ips(i,-2),th)+cps(i,-1)*ddr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*ddr(im,ir,ips(i,1),th)+cps(i,2)*ddr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns; do ir = -1, nx+2; do im = 1, mphi  
          ddr(im,ir,i,th) = aar(im,ir,i,th) * & 
            real(cps(i,-2)*bbr(im,ir,ips(i,-2),th)+cps(i,-1)*bbr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*bbr(im,ir,ips(i,1),th)+cps(i,2)*bbr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) - &  
            (ddr(im,ir-2,i,th)-8.0E0*ddr(im,ir-1,i,th)+8.0E0*ddr(im,ir+1,i,th)-ddr(im,ir+2,i,th))
        end do; end do; end do 

        ! third bracket 
        do i = 1, ns; do ir = 1, nx; do im = -1, mphi+2
          ddr(im,ir,i,th) = bbr(im,ir,i,th) * &
            (aar(im,ir-2,i,th)-8.0E0*aar(im,ir-1,i,th)+8.0E0*aar(im,ir+1,i,th)-aar(im,ir+2,i,th))
        end do; end do; end do 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi 
          ggr(im,ir,i,th) = ggr(im,ir,i,th) - &  
            real(cps(i,-2)*ddr(im,ir,ips(i,-2),th)+cps(i,-1)*ddr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*ddr(im,ir,ips(i,1),th)+cps(i,2)*ddr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns; do ir = -1, nx+2; do im = 1, mphi 
          ddr(im,ir,i,th) = bbr(im,ir,i,th) * & 
            real(cps(i,-2)*aar(im,ir,ips(i,-2),th)+cps(i,-1)*aar(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*aar(im,ir,ips(i,1),th)+cps(i,2)*aar(im,ir,ips(i,2),th))
        end do; end do; end do  
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) + &  
            (ddr(im,ir-2,i,th)-8.0E0*ddr(im,ir-1,i,th)+8.0E0*ddr(im,ir+1,i,th)-ddr(im,ir+2,i,th))
        end do; end do; end do  

        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ccr(im,ir,i,th) = ccr(im,ir,i,th) + rhostar * ggr(im,ir,i,th) * dtim *  &
                          & efun(ir,i,3,1) / (432.0E0 * mphi *dxgr *sgr_dist ) 
        end do; end do; end do 

        ! the first bracket 
        do i = 1, ns; do ir = 1, nx ; do im = 1, mphi   
          ggr(im,ir,i,th)= &
            real(cps(i,-2)*aar(im,ir,ips(i,-2),th)+cps(i,-1)*aar(im,ir,ips(i,-1),th)+                  &
            & cps(i,1)*aar(im,ir,ips(i,1),th)+cps(i,2)*aar(im,ir,ips(i,2),th)   ) *                &
            (bbr(im-2,ir,i,th)-8.0E0*bbr(im-1,ir,i,th)+8.0E0*bbr(im+1,ir,i,th)-bbr(im+2,ir,i,th))- &
            real(cps(i,-2)*bbr(im,ir,ips(i,-2),th)+cps(i,-1)*bbr(im,ir,ips(i,-1),th)+                  &
            & cps(i,1)*bbr(im,ir,ips(i,1),th)+cps(i,2)*bbr(im,ir,ips(i,2),th)   ) *                &
            (aar(im-2,ir,i,th)-8.0E0*aar(im-1,ir,i,th)+8.0E0*aar(im+1,ir,i,th)-aar(im+2,ir,i,th))
        end do; end do; end do    

        ! second bracket 
        do i = 1, ns; do ir = 1, nx;  do im = 1, mphi  
          ddr(im,ir,i,th) = aar(im,ir,i,th) * &
            (bbr(im-2,ir,i,th)-8.0E0*bbr(im-1,ir,i,th)+8.0E0*bbr(im+1,ir,i,th)-bbr(im+2,ir,i,th))
        end do; end do; end do; 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) + &  
            real(cps(i,-2)*ddr(im,ir,ips(i,-2),th)+cps(i,-1)*ddr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*ddr(im,ir,ips(i,1),th)+cps(i,2)*ddr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi 
          ddr(im,ir,i,th) = aar(im,ir,i,th) * & 
            real(cps(i,-2)*bbr(im,ir,ips(i,-2),th)+cps(i,-1)*bbr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*bbr(im,ir,ips(i,1),th)+cps(i,2)*bbr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) - &  
            (ddr(im-2,ir,i,th)-8.0E0*ddr(im-1,ir,i,th)+8.0E0*ddr(im+1,ir,i,th)-ddr(im+2,ir,i,th))
        end do; end do; end do 

        ! third bracket 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ddr(im,ir,i,th) = bbr(im,ir,i,th) * &
            (aar(im-2,ir,i,th)-8.0E0*aar(im-1,ir,i,th)+8.0E0*aar(im+1,ir,i,th)-aar(im+2,ir,i,th))
        end do; end do; end do 
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi 
          ggr(im,ir,i,th) = ggr(im,ir,i,th) - &  
            real(cps(i,-2)*ddr(im,ir,ips(i,-2),th)+cps(i,-1)*ddr(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*ddr(im,ir,ips(i,1),th)+cps(i,2)*ddr(im,ir,ips(i,2),th))
        end do; end do; end do 
        do i = 1, ns; do ir = -1, nx+2; do im = -1, mphi+2 
          ddr(im,ir,i,th) = bbr(im,ir,i,th) * & 
            real(cps(i,-2)*aar(im,ir,ips(i,-2),th)+cps(i,-1)*aar(im,ir,ips(i,-1),th) +  &
            & cps(i,1)*aar(im,ir,ips(i,1),th)+cps(i,2)*aar(im,ir,ips(i,2),th))
        end do; end do; end do  
        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ggr(im,ir,i,th) = ggr(im,ir,i,th) + &  
            (ddr(im-2,ir,i,th)-8.0E0*ddr(im-1,ir,i,th)+8.0E0*ddr(im+1,ir,i,th)-ddr(im+2,ir,i,th))
        end do; end do; end do  

        do i = 1, ns ; do ir = 1, nx; do im = 1, mphi
          ccr(im,ir,i,th) = ccr(im,ir,i,th) + rhostar * ggr(im,ir,i,th) * dtim *  &
                          & efun(ir,i,3,2) * lyinv / (432.0E0 *sgr_dist ) 
        end do; end do; end do 

      endif nlr


      ! do the back transform 
      do i = 1, ns; do ix = 1, nx
        call four1D_real(ccr(1:mphi,ix,i,th),aa(:,ix,i,th),FFT_FORWARD)
      end do; end do 
      
      if(erase_any_transfer_to) then
        do l = 1, size(no_transfer_to_modes)
          aa(no_transfer_to_modes(1,i),:, :, th) = 0
        end do
      end if
      
      ! get the size of 0,0 mode (should be tiny)
      ! reducing this over x could be expensive unless done outside loop
      ! not yet demonstrated that this is needed
      !dumc=sum(a(1,1:nx,th))/real(nx)
      !dumc=(0.,0.)

      ! add to rhs (may need to remove the 0,0 average)
      do ix = 1, nx; do imod = 1, nmod; do i = 1, ns 
        rhs(lincopy3(idx2,1,j,is)) = rhs(lincopy3(idx2,1,j,is)) - &
                                   & aa(imod,ix,i,th) 
        idx2 = idx2 + 1
      end do; end do; end do 

    end do vpar_loop
  end do loop3
  !!!$OMP end parallel do

end subroutine add_non_linear_terms_nonspec_2 



!------------------------------------------------------------------------------
!> This is a diagnostic routine. It appears here because it strongly relies 
!> on arrays that have been calculated in (and are local to) non-linear terms. 
!> This routine could of course be moved to radial diagnostic (needs some 
!> rewritting though. 
!>
!> The routine calculates the radial entropy profile as well as its flux. 
!------------------------------------------------------------------------------
subroutine entropy_radial(entr_rad) 

  use control,        only : nlapar
  use dist,           only : fmaxwl
  use dist,           only : fdisi, fdis_tmp2
  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp, n_x_grid, lsendrecv_x 
  use general,        only : gkw_abort
  use fft,            only : four1D_real, FFT_INVERSE
  use geom,           only : efun, bn, ints 
  use mode,           only : lyinv
  use components,     only : vthrat, signz, tmp 
  use velocitygrid,   only : vpgr, intvp, intmu
  use ompinterface,   only : ompget_thread_num
  use general,        only : gkw_abort
  use mpiinterface,   only : gather_array, mpiallreduce_sum
  use mpicomms,       only : COMM_DUMMY, COMM_X_EQ, COMM_X_NE

  real, intent(out) :: entr_rad(n_x_grid,2) 
  real, allocatable, save :: entr1(:,:), entr2(:,:)
  
  integer :: imod, ix, i, j, k, is, th, im, idx, idx2, iloop, n_fields
  integer :: ifields, ierr  
  real    :: ver 

  logical, save :: initialized = .false. 
  logical, save :: ALL_PROCS   = .true.  
  
  if (.not.initialized) then 
    allocate(entr1(nx,2),stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate entr1')
    allocate(entr2(nx,2),stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate entr2')
    initialized = .true. 
  endif 
  
  n_fields = 1
  if (nlapar) n_fields = 2

  ! outer loop over i, j, is 
  !do is = 1, nsp; do j = 1, nmu; do i = 1, ns  
 
   ! initialize to zero 
   entr_rad(:,:) = 0. 
   entr1(:,:) = 0. 
   entr2(:,:) = 0. 
 
  loop3: do iloop = 1, ns*nsp*nmu

    is = i3loop(1,iloop)
    j = i3loop(2,iloop)
    i = i3loop(3,iloop)

    ! Always 1 for non openmp runs
    th = ompget_thread_num() + 1
    idx=1

    ! fill the fields arrays
    do ifields = n_fields,1,-1 ! the linphi array returns iapar_ga first

      ! fill the array a with the potential 
      a(:,:,th) = (0.E0,0.E0) 

      if (lsendrecv_x) then 
        do ix = 1-xgp, nx+xgp; do imod = 1, nmod 
          a(imod,ix,th) = fdis_tmp2(linphi(idx,i,j,is))
          idx = idx + 1 
        end do; end do 
      else
        do ix = 1-xgp, nx+xgp; do imod = 1, nmod  
          a(imod,ix,th) = fdisi(linphi(idx,i,j,is))
          idx = idx + 1
        end do; end do
      endif 

     ! the boundary condition
      do imod = 1, nmod 
        a(imod,  -1,th) = mbd(-1,imod,i)*a(imod,ibd(-1,imod,i),th) 
        a(imod,   0,th) = mbd( 0,imod,i)*a(imod,ibd( 0,imod,i),th) 
        a(imod,nx+1,th) = mbd( 1,imod,i)*a(imod,ibd( 1,imod,i),th) 
        a(imod,nx+2,th) = mbd( 2,imod,i)*a(imod,ibd( 2,imod,i),th) 
      end do 
  
      ! do the inverse Fourier transform 
      ! FJC: Faster here to use fourcol -> fourcol_ra_ca ?
      ! which does all the rows of the array in 1D together
      do ix = 1-x_boundary_points, nx+x_boundary_points  
        call four1D_real(ar(1:mphi,ix,th),a(:,ix,th),FFT_INVERSE)
      end do
      
      ! Minimal time saving (~ 2%) and not working
      !call fourcol(ar(1:mphi,1-x_boundary_points:nx+x_boundary_points,th),&
      !     &a(:,1-x_boundary_points:nx+x_boundary_points,th),FFT_INVERSE)

      ! the periodic boundary conditions in the zeta direction 
      do ix = 1-x_boundary_points, nx+x_boundary_points 
        ar(    -1,ix,th) = ar(mphi-1,ix,th) 
        ar(     0,ix,th) = ar(  mphi,ix,th)
        ar(mphi+1,ix,th) = ar(     1,ix,th)
        ar(mphi+2,ix,th) = ar(     2,ix,th) 
      end do 

      ! store iapar_ga in er (this happens first time
      if (ifields == 2) er(:,:,th) =  ar(:,:,th)
      ! For electromagnetic runs the potential needs to be stored separately
      if (nlapar) then 
        if (ifields == 1) fr(:,:,th) = ar(:,:,th) 
      endif 

    end do ! ifields loop

    ! loop over the parallel velocity
    idx=1; idx2=1
    vpar_loop : do k = 1, nvpar
     
      ! copy out the distribution function 
      b(:,:,th) = (0.E0,0.E0)
      if (lsendrecv_x) then 
        do ix = 1-xgp, nx+xgp; do imod = 1, nmod 
          b(imod,ix,th) = fdis_tmp2(lindx3(idx,i,j,is))
          idx = idx + 1 
        end do; end do 
      else 
        do ix = 1-xgp, nx+xgp; do imod = 1, nmod
          !b(imod,ix,th) = fdis(indx(ifdis,imod,ix,i,j,k,is))
          b(imod,ix,th) = fdisi(lindx3(idx,i,j,is))
          idx = idx + 1
        end do; end do; 
      endif
      
      ! the boundary conditions 
      do imod = 1, nmod
        do ix = 1-x_boundary_points, 0
          b(imod,ix,th) = mbd(ix,imod,i)*b(imod,ibd(ix,imod,i),th) 
        end do
        do ix = 1, x_boundary_points
          b(imod,nx+ix,th) = mbd( 1,imod,i)*b(imod,ibd( 1,imod,i),th) 
        end do
      end do 

      ! do the inverse Fourier transform 
      do ix = 1-x_boundary_points, nx+x_boundary_points
        call four1D_real(br(1:mphi,ix,th),b(:,ix,th),FFT_INVERSE)
      end do 
      
      ! the periodic boundary conditions in the zeta direction 
      do ix = 1-x_boundary_points, nx+x_boundary_points 
        br(    -1,ix,th) = br(mphi-1,ix,th) 
        br(     0,ix,th) = br(  mphi,ix,th)
        br(mphi+1,ix,th) = br(     1,ix,th)
        br(mphi+2,ix,th) = br(     2,ix,th) 
      end do 

      ! put chi combined field into ar for electromagnetic runs 
      if (nlapar)  ar(:,:,th) = fr(:,:,th) - 2. * vpgr(i,j,k,is) * vthrat(is) * er(:,:,th)
 
      do ix = 1, nx 
        do im = 1, mphi;  
        
          ! the radial ExB velocity: \pd{\ga{phi}}{\zeta} \cdot \mathcal{E}^{\zeta\psi} 
          ver = (ar(im-2,ix,th)-8.E0*ar(im-1,ix,th)+8.E0*ar(im+1,ix,th)-ar(im+2,ix,th)) & 
              & * efun(ix,i,2,1) * lyinv / 12.E0 
              
          ! the entropy 
          entr1(ix,1) = entr1(ix,1) + ( br(im,ix,th)**2/(2.E0*fmaxwl(ix,i,j,k,is)) +    &
                      & signz(is)*br(im,ix,th)*ar(im,ix,th)/tmp(ix,is) )                &
                      & *intvp(i,j,k,is)*intmu(j)*bn(ix,i)*ints(i)
 
        ! the entropy flux 
          entr1(ix,2) = entr1(ix,2) + ( br(im,ix,th)**2/(2.E0*fmaxwl(ix,i,j,k,is)) +    &
                      & signz(is)*br(im,ix,th)*ar(im,ix,th)/tmp(ix,is) )                &
                      & *intvp(i,j,k,is)*intmu(j)*bn(ix,i)*ints(i)*ver 

        end do  
      end do 
 
    end do vpar_loop
  end do loop3
  
  ! global sum over all directions except the radial. 
  call mpiallreduce_sum(entr1,entr2,nx,2,COMM_X_EQ)
  
  ! then gather the array into one global array 
  call gather_array(entr_rad(1:n_x_grid,:),  n_x_grid,  2,           &
                  & entr2(1:nx, :),          nx,        2,           &
                  & COMM_X_NE, COMM_DUMMY, ALLGATHER = ALL_PROCS)
 

end subroutine entropy_radial




end module non_linear_terms



