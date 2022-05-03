!------------------------------------------------------------------------------
!> Sets up the GKW velocity grid. 
!> The nonuniform velocity grid is still under development.
!-------------------------------------------------------------------------------
module velocitygrid

  implicit none

  private
 
  ! publicly available routines 
  public :: velgrid_init, vgridboundary, vgridboundaryMom
  public :: elem_is_on_vpar_grid, get_vpar_stencil

  ! publicly available variables and parameters 
  public :: intmu,vpgr,mugr,dvp,dmu,intvp,iblow,ibhig
  public :: dvperp

  !The interval in the parallel velocity direction
  real, save :: dvp
  !The interval in mu
  real, save :: dmu
  !The interval in perpendicular velocity
  real, save :: dvperp 
  
  !> global rms parallel velocity
  real, public, save :: vpgr_rms
  !> global rms mu
  real, public, save :: mugr_rms

  !> the grid in parallel velocity space, vpgr(ns,nmu,nvpar,nsp) 
  real, allocatable, save :: vpgr(:,:,:,:)
  !> the grid in mu space,  mugr(nmu)
  real, allocatable, save :: mugr(:)

  !> the grid for velocity space integration, intmu(nmu) 
  real, allocatable, save :: intmu(:)
  !> the grid for velocity space integration, intvp(ns,nmu,nvpar,nsp)
  real, allocatable, save :: intvp(:,:,:,:)

  ! Arrays that determine at which position in the s-grid 
  ! a trapped particle bounces. (only used for vp_trap =1) 
  integer, allocatable, save :: iblow(:), ibhig(:)

  ! number of extra spaces needed in the parallel velocity direction
  integer, save :: ivpar_extra
  integer, save :: is_extra

  logical, save :: zero_grad = .false.

contains

subroutine velgrid_init

  call velgrid_allocate(1)
  call dist_grid_setup
  call nonuni_vel_grid_setup
  call velgrid_stats

end subroutine velgrid_init

!-----------------------------------------------------------------------------
!> This routine allocates all the arrays connected with the velocity grid
!> requires 1 call for the arrays used in parallel_plans,
!> then called with "1" to allocate and "-1" to deallocate
!-----------------------------------------------------------------------------

subroutine velgrid_allocate(i)

  use grid,    only : ns, nmu, nvpar, nsp
  use control, only : vp_trap, order_of_the_scheme
  use general, only : gkw_abort
  
  integer, intent(in) :: i 
  !  1 => allocate, 
  ! -1 => deallocate,

  ! the integer for the error message
  integer :: ierr

  ! initialize the error parameter 
  ierr= 0

  if (i .eq. 1) then

    ! allocate the theta grid array
    ! in some cases we need vpgr from regions we do not solve in
    if (order_of_the_scheme .eq. 'second_order') then
      ivpar_extra = 1
      is_extra = 1
      allocate(vpgr(0:ns+1,0:nmu+1,0:nvpar+1,nsp),stat=ierr)
      if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate vpgr')
      allocate(intvp(0:ns+1,0:nmu+1,0:nvpar+1,nsp),stat=ierr)
      if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate intvp')
    else if (order_of_the_scheme .eq. 'fourth_order') then 
      ivpar_extra = 2
      is_extra = 2
      allocate(vpgr(-1:ns+2,0:nmu+1,-1:nvpar+2,nsp),stat=ierr)
      if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate vpgr')
      allocate(intvp(-1:ns+2,0:nmu+1,-1:nvpar+2,nsp),stat=ierr)
      if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate intvp')
    else
      call gkw_abort('velgrid_allocate: ivpar_extra; bad case of scheme order')
      ivpar_extra = 0 
      is_extra = 0
      allocate(vpgr(ns,nmu,nvpar,nsp),stat=ierr)
      allocate(intvp(ns,nmu,nvpar,nsp),stat=ierr)
    endif
 
    intvp(:,:,:,:)=0.0
    vpgr(:,:,:,:) = 0.0
     
    ! allocate the velocity grid 
    allocate(mugr(0:nmu+1),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate mugr')
  
    ! allocate the array for velocity space integration
    ! mu-direction 
    allocate(intmu(0:nmu+1),stat=ierr)
    if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate intmu')
   
    if (vp_trap.eq.1) then 
       allocate(iblow(nvpar), stat = ierr) 
       if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate iblow')
       allocate(ibhig(nvpar), stat = ierr) 
       if (ierr.ne.0) call gkw_abort('dist_allocate: Could not allocate ibhig')
    endif
  
  else if (i .eq. -1) then

  ! deallocation

    if( allocated(mugr) )    deallocate(mugr)
    if( allocated(vpgr) )    deallocate(vpgr)
    if( allocated(intmu) )   deallocate(intmu)
    if( allocated(intvp) )   deallocate(intvp)

  else

    call gkw_abort('dist_allocate: called with i /= 1 or -1')

  endif

end subroutine velgrid_allocate

!-----------------------------------------------------------------------------
!> This routine calculates the grids in velocity space as well as the arrays
!> that are necessary for the integration over velocity space. The integral
!> over the distribution function should then be calculated as
!> B_N(i) sum_{jk} intmu(j) intvp(k) f(j,k)
!-----------------------------------------------------------------------------

subroutine dist_grid_setup 

  use grid,    only : n_mu_grid, lmu, ns, nmu, nsp 
  use grid,    only : n_vpar_grid, nvpar, gvpar, vpmax, mumax
  use control, only : uniform_mu_grid

  ! local parameters
  integer :: i, j, k, is
  real :: pi, vperp

  ! set the constant PI
  pi  = 4.E0*atan(1.E0)

  if (uniform_mu_grid) then 
    ! mu grid calculations  
    ! Calculate the mu grid values. 
    ! Calculate the help array for the mu integration 
    dvperp = sqrt(2.E0*mumax)/real(n_mu_grid)
    dmu = mumax /real(n_mu_grid)
    do j= 0,n_mu_grid+1
      if (lmu(j) .ge. 0 .and. lmu(j) .le. nmu+1) then
        mugr(lmu(j)) = (j-0.5)*dmu
        intmu(lmu(j)) = 2.E0*pi*dmu
      endif
    end do 
  else 
    dmu = mumax /real(n_mu_grid) !not used here...?
    dvperp = sqrt(2.E0*mumax)/real(n_mu_grid)
    do j = 0,n_mu_grid+1     
      vperp = (j-0.5)*dvperp 
      if (lmu(j) .ge. 0 .and. lmu(j) .le. nmu+1) then
        mugr(lmu(j)) = vperp**2 / 2.E0
        intmu(lmu(j)) =abs(pi*((vperp+0.5E0*dvperp)**2 - (vperp-0.5*dvperp)**2))
      endif 
    end do 
  endif 

  ! Caclulate the parallel velocity grid values 
  ! Calculate the help array for the parallel velocity integration 
  dvp = 2.0*vpmax / real(n_vpar_grid) 
  do is = 1, nsp
    do i = 1-is_extra, ns+is_extra 
      do j = 0, nmu+1
        do k = 1-ivpar_extra, nvpar+ivpar_extra
          vpgr(i,j,k,is) = - vpmax + (gvpar(k)-0.5E0)*dvp 
          intvp(i,j,k,is) = dvp 
        end do 
      end do
    end do 
  end do

end subroutine dist_grid_setup 

!-----------------------------------------------------------------------------
!> This routine sets up the velocity grid with arbitrary number of
!> points within the trapping condition determined by the value
!> of n_trapped. 
!> Must be called after parallelize_geom and parallelize_rotation
!-----------------------------------------------------------------------------

subroutine nonuni_vel_grid_setup
  use global,    only : r_tiny, gkw_a_equal_b_accuracy, int2char
  use control,   only : vp_trap
  use grid,      only : nvpar, ns, nmu, n_trapped, nperiod, N_s_grid, nsp
  use geom,      only : bn, bmin, bmax, geom_parallelized
  use rotation,  only : cfen, cfenh, cfenl, rotation_parallelized
  use specfun,   only : sf_erf 
  use constants, only : pi
  use gauss,     only : legendre_ek_compute
  use general,   only : gkw_abort

  real    :: b, vpmax, vpm, vpl,vpr, vpgr2, dvp2, blim, cfe
  real    :: ya, yb, val, xlow, xhig, xv, xhalf
  integer :: i, j, k, is, ierr, ndum, idum, m, kref, l, ix, nb 
  logical :: searching, verbose

  ! dummy array (low field side velocity) 
  real, allocatable :: vp0(:,:,:)
  real, allocatable :: vlftemp(:,:,:)

  ! dummy array (interval of the parallel velocity grid)
  real, allocatable :: dvp(:) 

  ! dummy array (selected points) 
  integer, allocatable :: ivps(:)

  !When number of bounce points is large, this is the array of 
  !values that will be kept.
  real, allocatable :: pkeep(:)

  ! Points and weights of the Gauss laguerre integration 
  real, allocatable :: x(:), w(:) 

  ! Here the value of ix is set and the calculation is (at the moment)
  ! done for one flux surface only 
  ix = 1 

  ! initial values
  vpl = 0.

  if (.not. geom_parallelized) then
    call gkw_abort('Vel grid setup: Geom not parallelized')
  end if
    
  if (.not. rotation_parallelized) then
    call gkw_abort('Vel grid setup: Rotation not parallelized')
  end if

  ! check whether the parallel velocity grid is to be set up with 
  ! the parallel velocity following the trapping condition 
  if (vp_trap.ne.1) return 

  ! Set the error parameter 
  ierr = 0 
  verbose = .false.

  ! set the number of bounce points
  if (mod(ns,2)==0) then
    nb = ns/2-1    ! even ns -> odd nb
  else 
    nb = (ns+1)/2  ! odd ns -> even nb
  end if 

  ! The set up in this routine assume that nperiod = 1. 
  if (nperiod.ne.1) call gkw_abort('With vp_trap = 1, nperiod must&
       & be equal to 1')

  ! nvpar must be even 
  if (mod(nvpar,2).ne.0) call gkw_abort('With vp_trap = 1 nvpar must&
       & be even')        
    
  if ((N_s_grid/2).lt.n_trapped) call gkw_abort('n_trapped is too big. &
       & The number of possible bounce points (ns/2-1) is smaller than  &
       & number of points in the trapped zone -> See manual')
  !The number of trapped points must be smaller than the number of cells 
  !in the parallel velocity direction
    
  if ((nvpar/2).lt.n_trapped) call gkw_abort('With vp_trap = 1 n_trapped must &
         & be smaller than the number of parallel grid points nvpar/2')
    
  if(((nvpar/2)-n_trapped).eq.1) call gkw_abort('You only have one point in   & 
         & the passing region of the velocity domain -> More are necessary.')

  ! Set the maximum parallel velocity (low field side) 
  vpmax = 3.E0 
  blim = 1.0E-13

  ! caclulate the low field side velocity for the trapped region 
  ! Only the positive values of v_parallel are caculated. These 
  ! are stored in vp0 
  allocate(vp0(nmu,nb+2,nsp),stat = ierr)
  if (ierr.ne.0) call gkw_abort('could not allocate vp0 in vel_grid_&
       &setup') 


  ! Determine the magnetic field at the bounce point of the 
  ! trapped particle, and the centrifugal energy
  ! (WARNING at present the latter only works for hydogenic plasma) 
  ! in which all speicies have same trapping condition
  ! Bounce points are defined exactly half way between points
  do j = 1, nmu; do is = 1, nsp
    do k = 1, nb+1  
      if (k.lt.nb+1) then 
        b = 0.5E0*(bn(ix,(ns+1)/2 + k) + bn(ix,(ns+1)/2 + k + 1)) 
        cfe = 0.5E0*(cfen((ns+1)/2 + k,is) + cfen((ns+1)/2 + k + 1,is)) 
      else 
        b = bmax(ix)
        cfe = cfenh(is) 
      endif
      ! the velocity on the low field side 
      vp0(j,k,is) = sqrt(2.E0*mugr(j)*(b - bmin(ix))+cfe-cfenl(is))
    end do
  end do; end do

  ! allocate the array for the storage of the v_parallel interval
  ! and the array that holds the selected points 
  allocate(dvp(nb+1), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Unable to allocate dvp in vel_grid&
       &setup') 
  allocate(ivps(nb+1), stat = ierr) 
  if (ierr.ne.0) call gkw_abort('Unable to allocate ivps in vel_grid&
       &setup')  
  allocate(vlftemp(nvpar/2,nmu,nsp), stat = ierr) !vlftemp is a temporary
  !array with the low field side velocities
  if (ierr.ne.0) call gkw_abort('Unable to allocate vlftemp in vel_grid&
       &setup')  
  allocate(pkeep(nb+1), stat = ierr)
  if (ierr.ne.0) call gkw_abort('Unable to allocate pkeep in vel_grid&
       &setup')  

  idum = 0
  do i=1, nb+1
    pkeep(i) = 1
    idum = idum +1
  end do
  write(*,*)'Initially there are', idum, 'points kept. N_trapped is',n_trapped
  ndum = nb+1-n_trapped
  write(*,*)ndum, 'Points must be neglected'       
  idum = 0

  !n_trapped points must now be selected for within the trapped region
  if (.true.) then 
  do while(idum.lt.ndum)
    do k = 2, nb
          
      searching = .true. 
      m = k - 1 
      do while(searching)
        if (m.eq.0) then 
          vpl = 0.E0 
          searching = .false.
        else 
          if (abs(pkeep(m)-1) < r_tiny) then 
            vpl = vp0(1,m,1) 
            searching = .false. 
          else 
            m = m - 1 
          endif
        endif
      end do
          
      searching = .true. 
      m = k + 1 
      do while (searching) 
        if (m.eq.nb+2) then 
          vpm = vp0(1,nb+2,1) 
          searching = .false. 
        else
          if (abs(pkeep(m) - 1) < r_tiny) then 
            vpm = vp0(1,m,1) 
            searching = .false. 
          else 
            m = m + 1 
          endif
        endif
      end do
          
      dvp(k) = vpm - vpl 
      if(verbose)then
        write(*,*)dvp(k),pkeep(k)
      end if
    end do
    !The point with the smallest interval is found and removed to even
    !out the resolution, and make sure n_trapped particles are kept.

    kref = 0 ; vpm = 1
    do j = 2, nb-1  
      if (abs(pkeep(j) - 1) < r_tiny) then 
        if (dvp(j).lt.vpm) then 
          kref = j 
          vpm = dvp(j) 
        endif
      endif
    end do
      
    ! eliminate the minimum point 
    pkeep(kref) = 0 

    !Count the number of points that have been removed 
    idum = idum + 1  
    if(verbose)then
      write(*,*)'idum = ',idum
    end if
  end do
  endif  ! End of section that is switched on / off   

  if(verbose)then
    do i=1,nb
      write(*,*)'Point along ns', i, 'Is it kept?' , pkeep(i)
    end do
  end if
   
  !Build the grid at the low field point -> This is the reference for
  !all other points around the Torus, vpar0
  do j = 1,nmu; do is = 1, nsp
       
    !The points inside the trapped region are determined by the velocity
    !of bounce points calculated earlier 
    idum = 0       
    do i=1,nb+1
      if(gkw_a_equal_b_accuracy(pkeep(i), 1.0)) then
        idum=idum+1
        vlftemp(idum,j,is) = vp0(j,i,is)
      end if
    end do

    if(idum.ne.n_trapped) call gkw_abort('Something problem')

    vlftemp(n_trapped+1,j,is) = 2.E0*vp0(j,nb+1,is)- vp0(j,nb,is) 
    !The first point outside the trapped region is equidistant from the
    !boundary as the first one inside
       
     dvp2 = ((vpmax - vlftemp(n_trapped+1,j,is))/(nvpar/2-n_trapped-1))

     !All other points outside this are evenly spaced
     write(*,*)'Velocity element outside trapped region = ',dvp2
     if(dvp2.lt.0)then
        call gkw_abort('Need more grid cells in the velocity grid outside &
             & the trapped region')
     end if
       
     if((nvpar/2).gt.n_trapped)then
       do i = n_trapped+2,nvpar/2
         vlftemp(i,j,is) = vlftemp(i-1,j,is) + dvp2
       end do
     end if
  end do; end do
    
  !Following commmented section prints the grid velocity grid points
  !to a file.
  open(18,file = 'velgridtemp') !The grid at the high field point
  do j = 1, nmu; do k =  1, nvpar/2 
    write(18,*) sqrt(mugr(j)),vlftemp(k,j,1),vp0(j,nb+1,1)   
  end do; end do

  deallocate(pkeep)

  !Section simultaneously calculates the velocity array along the field line
  !and symmetrises the array 
  do j = 1,nmu; do is = 1, nsp;  do i = 1, ns ; do k = 1,nvpar/2 
    vpgr2 = 2.E0*mugr(j)*(bmin(ix) - bn(1,i)) + cfenl(is) - cfen(i,is) + &
          & vlftemp(k,j,is)**2 
    if (vpgr2.gt.(0.E0)) then 
      vpgr(i,j,k+nvpar/2,is) = sqrt(vpgr2) 
      vpgr(i,j,nvpar/2-k+1,is) = -sqrt(vpgr2) 
      !vpgr(i,j,k+nvpar/2) = abs(vlftemp(k,j)) 
      !vpgr(i,j,nvpar/2-k+1) = - abs(vlftemp(k,j))
    else 
      vpgr(i,j,k+nvpar/2,is) = 0.E0
      vpgr(i,j,nvpar/2-k+1,is) = 0.E0 
    endif
  end do; end do; end do; end do
  deallocate(vlftemp)

  ! build the array for the velocity space integration 
  !It is implemented so that it ignores any cells set
  !to zero, finds the next non zero cell and calculates
  !interval from that.  Symmetrised on the fly as well.
  do i = 1, ns; do is = 1, nsp; do j = 1, nmu; do k = 1, nvpar/2 
    vpl = 0.E0
    vpr = 0.E0
    if (gkw_a_equal_b_accuracy(vpgr(i,j,k,is), 0.E0)) then
      intvp(i,j,k,is) = 0.E0 
      intvp(i,j,nvpar-k+1,is)=intvp(i,j,k,is)                 
    else 
      if (k.eq.1) then 
        vpl =  vpmax - abs(vpgr(i,j,k,is))
        kref = k+1
        searching = .true.
        do while(searching)
          if(abs(vpgr(i,j,kref,is)).lt.blim) then
            kref = kref+1
          else
            searching = .false.
          end if
        end do
        vpr = 0.5*abs((vpgr(i,j,kref,is)-(vpgr(i,j,1,is))))
        intvp(i,j,k,is) = vpr + vpl
        intvp(i,j,nvpar-k+1,is)=intvp(i,j,k,is)        
      else if (k.eq.(nvpar/2)) then                   
        vpr = abs(vpgr(i,j,nvpar/2,is))
        kref=k-1
        searching = .true.
        do while(searching)
          if(abs(vpgr(i,j,kref,is)).lt.blim)then
            kref=kref-1
          else
            searching= .false.
          end if
        end do
        vpl = 0.5*abs((vpgr(i,j,nvpar/2,is))-(vpgr(i,j,kref,is)))
        intvp(i,j,nvpar/2,is) = vpr+vpl
        intvp(i,j,nvpar-k+1,is)=intvp(i,j,k,is)    
      else
        kref = k+1
        searching = .true.
        do while(searching)
          if(abs(vpgr(i,j,kref,is)).lt.blim)then
            kref=kref+1
          else
            searching = .false.
          end if
        end do
        vpr = 0.5*abs((vpgr(i,j,kref,is)-(vpgr(i,j,k,is))))            
        kref=k-1
        searching = .true.
        do while(searching)
          if(abs(vpgr(i,j,kref,is)).lt.blim)then
            kref=kref-1
          else
            searching= .false.
          end if
        end do
        vpl = 0.5*abs((vpgr(i,j,k,is))- (vpgr(i,j,kref,is)))
        intvp(i,j,k,is) = vpr+vpl
        intvp(i,j,nvpar-k+1,is)=intvp(i,j,k,is)
      end if
    endif
  end do; end do; end do; end do

  !Following commmented section prints the grid velocity grid points
  !to a file.
  do is = 1, min(4,ns)
    open(18,file = 'velgridtemp'//int2char(is,1)) !The grid at the high field point
    do j = 1, nmu
      do k =  1, nvpar 
        write(18,*) sqrt(mugr(j)),vpgr(is,j,k,1),vp0(j,is,1)   
      end do
    end do
  end do

  ! switch off 
  if (.true.) then 

  do j = 1, nmu; do is = 1, nsp; do i = 1,ns 
    intvp(i,j,nvpar/2+n_trapped,is) = 0.5*intvp(i,j,nvpar/2+n_trapped,is) 
    intvp(i,j,nvpar/2-n_trapped+1,is) = intvp(i,j,nvpar/2+n_trapped,is) 
  end do; end do; end do 

  ! Now redo the passing particle gridding (GAUSSIAN INTEGRATION)
  m = nvpar/2 - n_trapped
  allocate(x(m), stat = ierr) 
  if (ierr /= 0) call gkw_abort('nonuni_vel_grid_setup: unable to &
                                 &allocate x') 
  allocate(w(m), stat = ierr) 
  if (ierr /= 0) call gkw_abort('nonuni_vel_grid_setup: unable to &
                                 &allocate w')

  do j = 1, nmu 

    do is = 1, nsp 

      ! transform the limits 
      ya = sf_erf(vp0(j,ns/2,is))  
      yb = sf_erf(vpmax) 

      ! calculate the points and weights of the integrator 
      call legendre_ek_compute (m, x, w) 

      do k = nvpar/2 + n_trapped + 1,  nvpar 
     
        m = k - nvpar/2 - n_trapped 
 
        ! determine the value of the velocity   
        val = ya + (yb-ya)*(1+x(m))/2 
        xlow = 0.
        xhig = vpmax 
        do l = 1, 52 
          xhalf = (xlow + xhig) / 2 
          !erf is fortran 2008. Does not compile on hector
          !xv    = erf(xhalf)
          xv    = sf_erf(xhalf)
          if (xv>val) then  
            xhig = xhalf 
          else 
            xlow = xhalf
          endif  
        end do 
        ! x now is the velocity (low field side) 
        x(m) = 0.5E0*(xlow + xhig) 
        w(m) = w(m)*exp(x(m)**2)*sqrt(pi)*(yb-ya)/4.E0 

        write(*,*)m,x(m),w(m)

      end do 

    end do 

    ! Then calculate the velocities and weights as a function 
    ! of the position along the field line 
    do i = 1, ns; do is = 1, nsp; do k = nvpar/2+n_trapped+1, nvpar 
     
      m = k - nvpar/2 - n_trapped 
      vpgr2 = 2.E0*mugr(j)*(bmin(ix) - bn(1,i)) + cfenl(is) - cfen(i,is) + x(m)**2 

      if (vpgr2 > 0) then 
        vpgr(i,j,k,is) = sqrt(vpgr2) 
        vpgr(i,j,nvpar-k+1,is) = - vpgr(i,j,k,is) 
        intvp(i,j,k,is) = w(m)*x(m)/vpgr(i,j,k,is)
        intvp(i,j,nvpar-k+1,is) = intvp(i,j,k,is)
      else  
        vpgr(i,j,k,is) = 0.
        vpgr(i,j,nvpar-k+1,is) = 0. 
        intvp(i,j,k,is) = 0.
        intvp(i,j,nvpar-k+1,is) = 0.
      endif 

    end do; end do; end do 

  end do ! loop over j (mu) 

  ! de-allocate the arrays of the points and weigths 
  deallocate(w) 
  deallocate(x)


  endif ! end of switch off 


  ! Finally fill the arrays that determine in which point in the 
  ! s-grid the trapped particles bounce 
  do k = 1, nvpar/2 - n_trapped  
    iblow(k) = 0 
    ibhig(k) = 0 
  end do
  do k = nvpar/2 + n_trapped + 1, nvpar
    iblow(k) = 0
    ibhig(k) = 0
  end do
  do k = nvpar/2 - n_trapped + 1, nvpar/2 + n_trapped 
    !write(*,*) k
    iblow(k) = 1 
    m = 1 
    do while (gkw_a_equal_b_accuracy(vpgr(m,1,k,1), 0.0) .and. (m.le.ns))
      m = m + 1 
      iblow(k) = m 
    end do
    if (m.gt.nb+1) call gkw_abort('vel_grid_setup: can not find iblow')
    ibhig(k) = ns 
    m = ns 
    do while(gkw_a_equal_b_accuracy(vpgr(m,1,k,1), 0.0) .and. (m.ge.1))
      m = m - 1 
      ibhig(k) = m 
    end do
    if (m.lt.nb+1) call gkw_abort('vel_grid_setup: can not find ibhig')
  end do

end subroutine nonuni_vel_grid_setup

!-----------------------------------------------------------------------------
!> calculate the rms values of the velocity grid
!-----------------------------------------------------------------------------

subroutine velgrid_stats

  use mpiinterface, only : mpiallreduce_sum, number_of_processors
  use grid,         only : ns,nmu,nvpar

  real, dimension(2) :: vpmu, buf
  
  !FJC: WARNING sort this out for multiple species!
  vpmu(1) = sum(vpgr(1:ns,1:nmu,1:nvpar,1)**2)/                      &
          & (1.*ns*nmu*nvpar*number_of_processors)
  vpmu(2) = sum(mugr(1:nmu)**2)/(1.*nmu*number_of_processors)

  call mpiallreduce_sum(vpmu,buf,2)

  vpgr_rms = sqrt(buf(1))
  mugr_rms = sqrt(buf(2))

end subroutine velgrid_stats

!-----------------------------------------------------------------------------
!>Function that deals with the boundary conditions for all terms
!>in the collision operator
!>ttype requires one of two options as an integer:
!> 1 - advective term
!> 2 - diffusive term
!---------------------------------------------------------------------------

subroutine vgridboundary(ttype,E,ingrid)

  use grid,       only : nvpar,nmu,n_vpar_grid,n_mu_grid
  use grid,       only : gvpar, gmu
  use structures, only : matrix_element
  use general,    only : gkw_abort

  !Short for term type.
  !If 1 then advective term -> Like friction boundary
  !If 2 then diffusive term -> Pitchangle etc
  integer, intent(in)                  :: ttype
  type (matrix_element), intent(inout) :: E
  logical, intent(out)                 :: ingrid

  integer :: kk,jj

  kk = E%kloc
  jj = E%jloc 

  if(ttype.ne.1)then
    if(ttype.ne.2)then
      call gkw_abort('ttype must be 1 or 2 in vgridboundary')
    end if
  end if

  if (gvpar(kk) .gt. n_vpar_grid) then
    ingrid = .false.
         
    if(zero_grad)then
      ingrid = .true.
      E%kloc = nvpar
    endif

    !Its a frictional term therefore the rest isnt needed
    if(ttype.eq.1)return

    !The corners of the grid will ALWAYS be set to zero.
    if(gmu(jj).lt.1)then
      ingrid=.false.
      return
    else if(gmu(jj).gt.n_mu_grid)then
      ingrid=.false.
      return
    else
      return
    end if
  else if (gvpar(kk) .lt. 1 ) then
    ingrid = .false.
   
    if(zero_grad)then
      E%kloc = 1
      ingrid = .true.
    endif
   
    if(ttype.eq.1)return
  
    if(gmu(jj).lt.1)then
      ingrid=.false.
      return
    else if(gmu(jj).gt.n_mu_grid)then
      ingrid=.false.
      return
    else
      return
    end if

  else if (gmu(jj).gt.n_mu_grid) then
    ingrid = .false.
   
    if(zero_grad)then
      E%jloc = nmu
      ingrid = .true.
    endif

    if(ttype.eq.1)return
   
    if(gvpar(kk).lt.1)then
      ingrid=.false.
      return
    else if(gvpar(kk).gt.n_vpar_grid)then
      ingrid=.false.
      return
    else
      return
    end if
 
  else if (gmu(jj).lt.1) then
   
    ingrid = .true.
    E%jloc = 1
   
    if(ttype.eq.1)return

    if(gvpar(kk).lt.1)then
      ingrid=.false.
      return
    else if(gvpar(kk).gt.n_vpar_grid)then
      ingrid=.false.
      return
    else
      return
    end if
  
  else
    ingrid = .true.
    return
  end if

  call gkw_abort('Error in vgridboundary')

end subroutine vgridboundary

!-----------------------------------------------------------------------------
!>Function that deals with the boundary conditions for all terms
!>in the collision operator
!>ttype requires one of two options as an integer:
!> 1 - advective term
!> 2 - diffusive term
!-----------------------------------------------------------------------------

subroutine vgridboundaryMom(ttype,jj,kk,j0,k0,ingrid)

  use grid,      only : nvpar,nmu, n_vpar_grid, gvpar, n_mu_grid, gmu
  use general,   only : gkw_abort

  !Short for term type.
  !If 1 then advective term -> Like friction boundary
  !If 2 then diffusive term -> Pitchangle etc
  integer, intent(in)    :: ttype
  integer, intent(inout) :: jj,kk
  integer, intent(in)    :: j0,k0
  !The prefactors for the mu boundary conditions.
  logical, intent(out)   :: ingrid

  if(ttype.ne.1)then
    if(ttype.ne.2)then
      call gkw_abort('ttype must be 1 or 2 in vgridboundary')
    end if
  end if

  !Extermal boundaries (i.e. ones at the edge of the computational
  !domain
  if (gvpar(k0).eq.(n_vpar_grid+1)) then
    ingrid = .false.
    return
  
  else if (gvpar(k0).eq.0)then 
    ingrid = .false.
    return

  else if (gmu(j0).eq.(n_mu_grid+1)) then
    ingrid = .false.
    return
  
  else if (gmu(j0).eq.0)then 
    ingrid = .false.
    return
    
  !Internal boudaries (i.e. processor boundaries)
  else if (kk.gt.nvpar) then
    ingrid = .false.
         
    if(zero_grad.and.(gvpar(kk).gt.n_vpar_grid))then
      ingrid = .true.
      kk=nvpar
    endif

    !Its a frictional term therefore the rest isnt needed
    if(ttype.eq.1)return

    !The corners of the grid will ALWAYS be set to zero.
    if(gmu(jj).lt.1)then
      ingrid=.false.
      return
    else if(gmu(jj).gt.n_mu_grid)then
      ingrid=.false.
      return
    else
      return
    end if
  else if (kk .lt. 1) then
    ingrid = .false.
   
    if(zero_grad.and.(gvpar(kk).lt.1))then
      kk=1  
      ingrid = .true.
    endif
   
    if(ttype.eq.1)return

    if(gmu(jj).lt.1)then
      ingrid=.false.
      return
    else if(gmu(jj).gt.n_mu_grid)then
      ingrid=.false.
      return
    else
      return
    end if

  else if (jj .gt. nmu) then
    ingrid = .false.
   
    if(zero_grad.and.(gmu(jj).gt.n_mu_grid))then
      jj = nmu
      ingrid = .true.
    endif

    if(ttype.eq.1)return
   
    if(gvpar(kk).lt.1)then
      ingrid=.false.
      return
    else if(gvpar(kk).gt.n_vpar_grid)then
      ingrid=.false.
      return
    else
      return
    end if
  else if (jj .lt. 1) then
    ingrid = .false.   

    if(gmu(jj).lt.1)then
      ingrid = .true.
      jj=1
    endif
   
    if(ttype.eq.1)return

    if(gvpar(kk).lt.1)then
      ingrid=.false.
      return
    else if(gvpar(kk).gt.n_vpar_grid)then
      ingrid=.false.
      return
    else
      return
    end if
  else
    ingrid = .true.
    return
  end if
  call gkw_abort('Error in vgridboundaryMom')

end subroutine vgridboundaryMom

!-----------------------------------------------------------------------------
!> This function decides if a central or a more or less one-sided stencil
!> is used.
!-----------------------------------------------------------------------------
function get_vpar_stencil(E) result(ist)
  use control, only : lflapv
  use structures, only : matrix_element
  use grid, only : gvpar, n_vpar_grid
  type (matrix_element), intent(in) :: E
  !> true is returned if the element is indeed on the parallel velocity grid.
  integer :: ist

  if(lflapv)then
    if(gvpar(E%kloc) == 1)then
      ist = -2
    else if(gvpar(E%kloc) == 2)then
      ist = -1
    else if(gvpar(E%kloc) == n_vpar_grid)then
      ist = 2
    else if(gvpar(E%kloc) == n_vpar_grid-1)then
      ist = 1
    else
      ist = 0
    end if

  else
    ist = 0
  end if

end function get_vpar_stencil

!-----------------------------------------------------------------------------
!> This function determines if the column index of the given matrix element is
!> within the parallel velocity grid.
!-----------------------------------------------------------------------------
function elem_is_on_vpar_grid(E)
  use general, only : gkw_abort
  use control, only : vp_trap
  use grid, only : nvpar, gvpar, n_vpar_grid
  use structures, only : matrix_element
  type (matrix_element), intent(in) :: E
  logical :: elem_is_on_vpar_grid

  ! the default - the other cases below must be .true.
  elem_is_on_vpar_grid = .false.
  ! always in the vpar grid in this case
  if (E%kloc >= -1 .and. E%kloc <= (nvpar+2)) elem_is_on_vpar_grid = .true.
  !It goes from -1 to +2 because there are two ghost points
  !because of the 4th order scheme.  2nd order there should only
  !be one point either side....dont know if this is needed though.

  ! in parallel_vpar case (with vp_trap = 0), we are in the grid
  ! for  k < 1, unless it is the first processor in vpar and for
  ! k > nvpar, unless we are the last processor in vpar (i.e. on
  ! the boundary) Note that the minimum and maximum allowed k
  ! below depend on the order of the scheme and should really be
  ! set elsewhere.
  if (vp_trap == 0) then

    !Only the points off the global grid are out of elem_is_on_vpar_grid=.false.
    if(gvpar(E%kloc) < 1) elem_is_on_vpar_grid = .false.
    if(gvpar(E%kloc) > n_vpar_grid) elem_is_on_vpar_grid = .false.

    if (E%kloc < 1-2 .or. E%kloc > nvpar + 2) then
      call gkw_abort('elem_is_on_vpar_grid: something is wrong')
    end if
  end if
end function elem_is_on_vpar_grid

end module velocitygrid
