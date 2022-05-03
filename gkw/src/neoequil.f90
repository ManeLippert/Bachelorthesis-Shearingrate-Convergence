!----------------------------------------------------------------------
!> Contains the interface with NEO and the calculation of neoclassical 
!> equilibrium components
!---------------------------------------------------------------------

module neoequil

  use global,       only : lenswitch

  implicit none

  private

  public :: neoequil_init, neoequil_check_params, rrefoa
  public :: gradneof, neof, mid_r_point, neophi
  public :: gradneophi, neo_eps

  ! public for diagnostic output
  public :: neo_coeffs_inter

  !> a list of integers, to be specified in the input file;
  !> The leading elements must all be > 0, and the rest -1. It is not allowed
  !> to have a -1 in between somewhere.
  integer, save, public :: neo_equil_parse_sp_seq(32) = -1

  !> The neoclassical distribution function,  neof(n_species,n_e,n_xi,n_theta) 
  real, allocatable, save :: neof(:,:,:,:,:)
  real, allocatable, save :: neophi(:,:)
  real, allocatable, save :: neophi_G(:,:)
  real, allocatable, save :: neophi_inter(:,:)
  real, allocatable, save :: gradneof(:,:,:,:)
  real, allocatable, save :: gradneophi(:)
  real, allocatable, save :: neo_coeffs(:,:,:,:,:)
  real, allocatable, save :: neo_coeffs_inter(:,:,:,:,:)
  real, allocatable, save :: neo_equilparam(:,:)
  real, allocatable, save :: neo_normparam(:,:)

  !> The theta points used in NEO
  real, allocatable, save :: neo_theta(:), neo_eps(:)

  !> The gridsizes of the distribution function from NEO
  integer, save :: neo_n_theta, neo_n_species, neo_n_energy, neo_n_xi
  integer, save :: neo_n_radial, mid_r_point
  real   , save :: rrefoa

  !> is the index of the first -1 element in neo_equil_parse_sp_seq,
  !> is the number of GKW species that get NEO equilibria
  !> ("neo_n_species incl. duplicates"). this is <= number_of_species,
  !> but can be larger or smaller than neo_n_species.
  integer, save, public :: neo_n_species_incl_dup

  !> Length of chararacter string of filename 
  integer, parameter :: lenfile = 180

  !> the first neo_nsp species on the local process have NEO
  !> equilibria. Of course, neo_nsp <= nsp always. Should be zero on processes
  !> whose species do not have NEO equilibria.
  integer, save, public :: neo_nsp
contains

  !------------------------------------------------------------------------------
  !>
  !------------------------------------------------------------------------------
  subroutine neoequil_check_params
  
  use control, only : flux_tube
  use general, only : gkw_abort
  use grid,    only : nperiod
  
  logical :: neofile_exists  
  character(len = lenfile) :: neofile
 
  !Check that these, essential, files are present so that they can be read.
  neofile = "out.neo.f"
  inquire(file=neofile,EXIST=neofile_exists)
  if (.not. neofile_exists) call gkw_abort('out.neo.f: not found: '//neofile)

  neofile = "out.neo.grid"
  inquire(file=neofile,EXIST=neofile_exists)
  if (.not. neofile_exists) call gkw_abort('out.neo.grid: not found: '//neofile)
  if(.not.flux_tube) call gkw_abort('Neoclassical correction only implemented &
     & for flux tube geometry')
  if(nperiod > 1) call gkw_abort('Neoclassical correction still requires nperiod = 1')

  neofile = "out.neo.phi"
  inquire(file=neofile,EXIST=neofile_exists)
  if (.not. neofile_exists) call gkw_abort('out.neo.phi: not found: '//neofile)

  neofile = "out.neo.equil"
  inquire(file=neofile,EXIST=neofile_exists)
  if (.not. neofile_exists) call gkw_abort('out.neo.equil: not found: '//neofile)

  neofile = "out.neo.expnorm"
  inquire(file=neofile,EXIST=neofile_exists)
  if (.not. neofile_exists) call gkw_abort('out.neo.expnorm: not found: '//neofile)

  end subroutine neoequil_check_params

  !------------------------------------------------------------------------------
  !>
  !------------------------------------------------------------------------------
  subroutine neoequil_allocate

  use grid,    only : ns
  use grid,    only : nmu, nvpar
  use grid,    only : nsg_pt
  use general, only : gkw_abort
  use control, only : order_of_the_scheme  
  use mpiinterface, only : root_processor

  integer :: ierr
  
  ! intialize the error parameter
  ierr=0
  
  !First column is only 1:5 as this is the number of radial grid points needed
  !to calculate the radial gradients
  if (order_of_the_scheme .eq.'fourth_order') then 
    if (neo_nsp > 0) then
      allocate(neof(1:5,1:neo_nsp,-1:ns+2,1:nmu,-1:nvpar+2),stat=ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')
      neof = 0.0
    end if

    allocate(neophi(1:5,-1:ns+2),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate neophi in neoequil')
  else if (order_of_the_scheme .eq.'second_order') then
    if (neo_nsp > 0) then
      allocate(neof(1:5,1:neo_nsp,0:ns+1,1:nmu,0:nvpar+1),stat=ierr)
      if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')
      neof = 0.0
    end if

    allocate(neophi(1:5,0:ns+1),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate neophi in neoequil')
    neophi(:,:) = 0.E0
  else
    call gkw_abort('velgrid_allocate: ivpar_extra; bad case of scheme order')
  endif

  if (neo_nsp > 0) then
    allocate(gradneof(1:neo_nsp,1:ns,1:nmu,1:nvpar),stat=ierr)
    if (ierr /= 0) call gkw_abort('Could not allocate gradneof in neoequil')

  end if
  allocate(gradneophi(1:ns),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate gradneophi in neoequil')
  gradneophi = 0.0
  
  allocate(neophi_G(1:5,1:neo_n_theta+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neophi in neoequil')
  neophi(:,:) = 0.E0
 
  allocate(neophi_inter(1:5,-1:nsg_pt+2),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neophi_inter in neoequil')
  neophi_inter(:,:) = 0.E0

  allocate(neo_coeffs(1:neo_n_radial,1:neo_n_species,  &
     & 1:neo_n_energy+1,1:neo_n_xi+1,1:neo_n_theta+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')
  neo_coeffs(:,:,:,:,:) = 0.E0

  allocate(neo_coeffs_inter(1:neo_n_radial,1:neo_n_species,-1:nsg_pt+2,  &
     & 1:neo_n_energy+1,1:neo_n_xi+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')
  neo_coeffs(:,:,:,:,:) = 0.E0
  
  allocate(neo_equilparam(1:neo_n_radial,1:(7+5*neo_n_species)),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_equilparam in neoequil')
  neo_equilparam(:,:) = 0.E0

  allocate(neo_normparam(1:neo_n_radial,1:(7*neo_n_species)),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_normparam in neoequil')
  neo_normparam(:,:) = 0.E0

  if(root_processor)write(*,*)'Neo arrays allocated'

end subroutine neoequil_allocate

!------------------------------------------------------------------------------
!> deallocate arrays which are only needed for initialisation
!------------------------------------------------------------------------------
subroutine neoequil_deallocate_tmps

  if(allocated(neo_coeffs)) deallocate(neo_coeffs)
  if(allocated(neo_coeffs_inter)) deallocate(neo_coeffs_inter)
  !if(allocated(neo_eps)) deallocate(neo_eps)
  if(allocated(neophi_G)) deallocate(neophi_G)
  if(allocated(neophi_inter)) deallocate(neophi_inter)
  if(allocated(neo_theta)) deallocate(neo_theta)
  if(allocated(neo_equilparam)) deallocate(neo_equilparam)
  if(allocated(neo_normparam)) deallocate(neo_normparam)

end subroutine neoequil_deallocate_tmps

subroutine read_gkw_neo

  use io,           only : get_free_file_unit

  character(len = lenfile) :: filename2
  integer :: file_unit
  logical :: lread_input

  filename2 = "FDS.dat"
   
  inquire(FILE=filename2,EXIST=lread_input)
  if (lread_input) then
     call get_free_file_unit(file_unit)
     open(file_unit,file=filename2,FORM='formatted',STATUS='unknown')
     close(file_unit)
  endif

end subroutine read_gkw_neo

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine read_neo_gridsizes

  use mpiinterface, only : root_processor, mpibcast
!  use mpiinterface, only : processor_number
  use io,           only : get_free_file_unit
  use general,      only : gkw_abort, gkw_warn
  use grid,         only : number_of_species, n_s_grid, nsp, lsp, gsp
  use components,   only : signz_G
  use global, only : int2char, r_tiny

  character(len = lenfile) :: neofile2
  integer :: j, file_unit, ierr
  integer :: i


  !Check that these, essential, files are present so that they can be read.
  neofile2 = "out.neo.grid"
  if(root_processor) write(*,*) 'Reading from NEO output files'
  
  ! Open ASCII file
  if(root_processor)then
    call get_free_file_unit(file_unit)
    open(file_unit, status='unknown',action='read', file=neofile2)  
    !Read the number of species
    read(file_unit,*) neo_n_species

    if(neo_n_species > number_of_species)then
      write(*,*) 'Neo ', neo_n_species, ' GKW ', number_of_species 
      call gkw_warn('The number of species from NEO is greater than that &
         & of GKW.')
    endif

    ! count, how many elements in the parse sequence are > 0
    neo_n_species_incl_dup = 0
    ! could loop only until number_of_species, but to spot errors, check more:
    do neo_n_species_incl_dup = 1, size(neo_equil_parse_sp_seq)
      if(neo_equil_parse_sp_seq(neo_n_species_incl_dup) < 0) then
        exit
      end if
    end do
    neo_n_species_incl_dup = neo_n_species_incl_dup - 1
    ! if no elements are specified this means that the default sequence
    ! shall be used: no species get the same background
    if(neo_n_species_incl_dup == 0) then
      call gkw_warn('Will use the default sequence to parse neoequil &
         & species: 1,2,3,...')
      neo_n_species_incl_dup = neo_n_species
      do i = 1, neo_n_species
        neo_equil_parse_sp_seq(i) = i
      end do
      write (*,*) neo_equil_parse_sp_seq(1:neo_n_species)
    end if
    if(neo_n_species_incl_dup < neo_n_species) then
      call gkw_warn('Only '//int2char(neo_n_species_incl_dup)//' of ' &
         & //int2char(neo_n_species)// ' NEO species are used.')
    end if
    if(neo_n_species_incl_dup > number_of_species) then
      call gkw_abort('There are more elements in neo_equil_parse_sp_seq &
         & than GKW species are defined.')
    end if
    if(maxval(neo_equil_parse_sp_seq) > neo_n_species) then
      call gkw_abort('The element '//int2char(maxval(neo_equil_parse_sp_seq))//&
         & ' appears in neo_equil_parse_sp_seq, but there are only ' &
         & //int2char(neo_n_species)//' species in the NEO data.')
    end if

  
    !Read the number of energy points
    read(file_unit,*) neo_n_energy
 
    !Read the number of pitch-angle points
    read(file_unit,*) neo_n_xi

    !Read the number of poloidal grid points
    read(file_unit,*) neo_n_theta
      
    if(neo_n_theta.ne.n_s_grid)then
      write(*,*) 'Neo: number of theta points', neo_n_theta 
      call gkw_warn('The number of s points is not equal')
    endif
  endif

  call mpibcast(neo_equil_parse_sp_seq,size(neo_equil_parse_sp_seq))
  call mpibcast(neo_n_species_incl_dup,1)
  call mpibcast(neo_n_species,1)
  call mpibcast(neo_n_energy, 1)
  call mpibcast(neo_n_xi,1)
  call mpibcast(neo_n_theta,1)

  ! the first neo_nsp local species have NEO backgrounds
  neo_nsp = lsp(min(gsp(nsp), neo_n_species_incl_dup))
  ! write (*,*) processor_number, " has neo_nsp", neo_nsp, nsp, &
  !   & neo_n_species_incl_dup

  !Values of theta !The plus one is to add a +pi point (-pi is there)
  !On all processors
  allocate(neo_theta(1:neo_n_theta+1),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')

  if(root_processor)then
    do j=1,neo_n_theta
      read(file_unit,*) neo_theta(j)
    end do
    neo_theta(neo_n_theta+1)=-neo_theta(1) 
  endif

  call mpibcast(neo_theta,neo_n_theta+1)

  if(root_processor)then
    !Number of radial grid points when NEO is run in profile mode 
    read(file_unit,*) neo_n_radial
  endif

  call mpibcast(neo_n_radial,1)

  if(neo_n_radial.ne.5)then
    if(root_processor)write(*,*) 'Neo: number of radial points', neo_n_radial 
    call gkw_warn('The number of radial points from NEO should be 5')
    if(neo_n_radial.eq.1)call gkw_warn('1 Radial grid point.  Radial gradients ignored')
  endif

  allocate(neo_eps(1:neo_n_radial),stat=ierr)
  if (ierr /= 0) call gkw_abort('Could not allocate neo_theta in neoequil')

  if(root_processor)then
    do j=1,neo_n_radial
      read(file_unit,*) neo_eps(j)
    end do
    close(file_unit)
  endif
  call mpibcast(neo_eps,neo_n_radial)

  if(neo_n_radial.eq.1) then
    mid_r_point = 1
  else
    mid_r_point = (neo_n_radial+1)/2
  end if


  if(root_processor)then
    write(*,*) 'Neo rad ,', 'Neo species ,', 'N energy ,', 'N xi ,', 'N theta'  
    write(*,*) neo_n_radial, neo_n_species, neo_n_energy, neo_n_xi, neo_n_theta
    if(abs(signz_G(maxloc(neo_equil_parse_sp_seq,1)) - (-1)) > r_tiny) then
      ! this test supposedly does not work correctly if several species
      ! have the NEO electron background
      call gkw_abort('The NEO electrons (last species) are not associated to &
       & the GKW electrons.')
    endif
  endif

end subroutine read_neo_gridsizes

!------------------------------------------------------------------------------
!> Read coefficients and normalisation parameters.
!------------------------------------------------------------------------------
subroutine read_neo_f
  use mpiinterface, only : root_processor, mpibcast
  use io,           only : get_free_file_unit

  character(len = lenfile) :: neofile
  integer :: j, file_unit, nelements
  integer :: is, i, k, ix, ncol

  neofile = "out.neo.f"

  nelements = neo_n_radial*neo_n_species*(neo_n_energy+1)*(neo_n_xi+1)*neo_n_theta
  
  if(root_processor) then
    call get_free_file_unit(file_unit)
    open(file_unit, status='unknown',action='read', file=neofile)
    do ix = 1,neo_n_radial
      do is = 1,neo_n_species
        do i = 1,neo_n_energy+1
          do j = 1,neo_n_xi+1
            do k = 1,neo_n_theta
              !Put the coefficients in a structure that is a little more managable.
              read(file_unit,*) neo_coeffs(ix,is,i,j,k)
            enddo
          enddo
        enddo
      enddo
    enddo
    close(file_unit)
  endif
  call mpibcast(neo_coeffs, nelements)

  !Read the equil file, to find the normalising quantities
  if(root_processor)then
    call get_free_file_unit(file_unit)
    neofile = "out.neo.equil"
    open(file_unit, status='unknown',action='read', file=neofile)  
    ncol = 7 + 5*neo_n_species
 
    do k = 1,neo_n_radial
       read(file_unit,*) (neo_equilparam(k,i),i=1,ncol)
    enddo  
  endif

  call mpibcast(neo_equilparam,  (7+5*neo_n_species)*neo_n_radial)

  !The ratio of Rref/a needed for transforming between NEO and GKW units
  rrefoa = neo_equilparam(1,5)

  !Read the equil file, to find the normalising quantities
  if(root_processor)then
    call get_free_file_unit(file_unit)
    neofile = "out.neo.expnorm"
    open(file_unit, status='unknown',action='read', file=neofile)  
    ncol = 7
 
    do k = 1,neo_n_radial
       read(file_unit,*) (neo_normparam(k,i),i=1,ncol)
    enddo  
  endif

  call mpibcast(neo_normparam,7*neo_n_radial)

end subroutine read_neo_f

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine read_neo_phi

  use mpiinterface, only : root_processor, mpibcast
  use io,           only : get_free_file_unit

  character(len = lenfile) :: neofile
  integer :: file_unit, ncol, nrow = 0
  integer :: i, j, k

  !Check that these, essential, files are present so that they can be read.
  neofile = "out.neo.phi"
    
  ! Open ASCII file
  if(root_processor)then
    call get_free_file_unit(file_unit)
    open(file_unit, status='unknown',action='read', file=neofile)  
    ncol = 3
    !This whole function assumes that NEO phi ouput is column ordered data
    !with 3 columns and n rows
    
    if(modulo(neo_n_theta,ncol).eq.0)then
      nrow = neo_n_theta/ncol    
    else if(modulo(neo_n_theta,ncol).eq.1)then
      nrow = (neo_n_theta+2)/ncol    
    else if(modulo(neo_n_theta,ncol).eq.2)then    
      nrow = (neo_n_theta+1)/ncol
    end if
    
    do k = 1,neo_n_radial
      ncol = 3
      do j = 1,nrow
        if(j.eq.nrow)ncol=modulo(neo_n_theta,ncol)
        read(file_unit,*) (neophi_G(k,3*(j-1)+i),i=1,ncol)
      enddo
      neophi_G(k,neo_n_theta+1)=neophi_G(k,1)
    enddo

  endif

  call mpibcast(neophi_G, neo_n_radial*(neo_n_theta+1) )  
  
end subroutine read_neo_phi

!------------------------------------------------------------------------------
!> Reconstructs the distribution function and transforms into
!> GKW velocity space
!------------------------------------------------------------------------------
subroutine reconstruct_f
  
  use grid,      only : nvpar, ns, nmu, n_s_grid
  use grid,      only : nmu
  use grid,      only : gs, gsp
  use geom,      only : bn, dBdpsi
  use velocitygrid,   only : vpgr, mugr
  use functions, only : legendre, leguerre
  use components, only : signz

  real    :: zeta, ener, ener_scal, kcoef, dx
  integer :: ne,nxi,is,i,j,k,ix,ii,iii
  real    :: bloc
  integer :: neo_is

  do ix=1,neo_n_radial
    ! do something for every local species for which there is a NEO
    ! equilibrium
    do is = 1, neo_nsp
      do i = -1,ns+2

        neo_is = neo_equil_parse_sp_seq(gsp(is))

        !The reconstructed distribution function (and its radial
        !gradient is stored locally on each processor, so this 
        !global s point is needed
        ii = gs(i)
        neophi(ix,i) = neophi_inter(ix,ii) 
        neophi(ix,i) = neophi(ix,i)*neo_normparam(ix,5)/neo_normparam(3,5)
        
        do j = 1,nmu
          do k = -1,nvpar+2
          
            !The corresponding pitch angle and energy value for 
            !specific GKW velocity grid point
            dx = neo_eps(ix)-neo_eps(3)
            !Magnetic field is Taylor expanded and linearly interpolated
            !To other surfaces
            if(ii.eq.0)then
               iii=n_s_grid
            elseif(ii.eq.-1)then
               iii=n_s_grid-1
            elseif(ii.eq.n_s_grid+1)then
               iii=1
            elseif(ii.eq.n_s_grid+2)then
               iii=2
            else
               iii=ii
            endif
            bloc = bn(1,i) + dBdpsi(1,iii)*dx/rrefoa

            !Pitch angle and energy coordinates
            zeta =  vpgr(1,j,k,is)/sqrt(vpgr(1,j,k,is)**2 + 2.E0*mugr(j)*bloc)
            ener =  (vpgr(1,j,k,is)**2 + 2.E0*mugr(j)*bloc)

            !Rescale the reference thermal velocity
            ener = ener*(neo_normparam(3,6)/neo_normparam(ix,6))**2
           
            !Rescales the energy so that the argument of the
            !Chebyshev polynomial is 2sqrt(en/en_max) - 1
            !ener_scal = 2*sqrt(ener) - 1 !!Not used any more.

            !Here we sum over the polynomials to reconstruct f in 
            !GKWs velocity space coordinates   
            neof(ix,is,i,j,k)=0.E0
            do ne = 1,neo_n_energy+1
              do nxi = 1,neo_n_xi+1
                if(nxi.eq.1)then
                  kcoef = 0.5E0
                  ener_scal = 1.0E0
                else
                  kcoef = 1.5E0
                  ener_scal = sqrt(ener)
                endif         
                !This expansion is used in the NEO 2009 PPCF paper
                !neof(ix,is,i,j,k) = neof(ix,is,i,j,k) + neo_coeffs(ix,is,ne,nxi,i)* & 
                 !& legendre(nxi-1,zeta)*cheby(ne-1,ener_scal)
                !This one is used in the 2012 PPCF
                !write(*,*) ne-1, nxi-1, i, j, k, neo_coeffs_inter(ix,is,ii,ne,nxi)
                neof(ix,is,i,j,k) = neof(ix,is,i,j,k) + &
                   & neo_coeffs_inter(ix,neo_is,ii,ne,nxi)* & 
                 & legendre(nxi-1,zeta)*ener_scal*leguerre(ne-1,ener,kcoef)
              enddo
            enddo 
            
            !Remove the adiabatic component
            !tspec = neo_equilparam(ix,7+is*2)           
            neof(ix,is,i,j,k) = neof(ix,is,i,j,k) &
                & - signz(is)*neophi(ix,i)
                       
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine reconstruct_f

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine neof_renormalise

  use grid,      only : nvpar, ns, nmu
  use grid,      only : nmu

  integer :: is,i,j,k,ix

  !The normalisations here are due to the fact that we use 5 radial neo calculations
  !and the normalising quantaties vary radially.  It is assumed that the middle one has
  !the same values as the GKW normalisations.
  do ix=1,neo_n_radial
    do i = -1,ns+2
      !Renormalise to the GKW rhostar 
      neophi(ix,i) = sqrt(2.E0)*neophi(ix,i)/rrefoa
      do j = 1,nmu
        do k = -1,nvpar+2
          do is = 1, neo_nsp
            
            !Renormalise to the GKW rhostar 
            neof(ix,is,i,j,k) = sqrt(2.E0)*neof(ix,is,i,j,k)/rrefoa

          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine neof_renormalise

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine cal_gradients

  use grid,      only : nvpar, ns, nmu
  use grid,      only : nmu

  integer :: is,i,j,k
  real :: dum,dx

  !NEO is run in profile mode.  5 radial points are taken with the outer
  !4 used to calculate the radial gradient of the distribution function
  !and the middle one used for the value.
  if(neo_n_radial.eq.1)then
    gradneof(:,:,:,:)=0.E0
    return
  endif

  dx = neo_eps(2)-neo_eps(1)
  
  do i = 1,ns

    !Second order finite difference radially
    dum = 8*(neophi(4,i)-neophi(2,i))
    dum = dum + (-neophi(5,i)+neophi(1,i))
    gradneophi(i) = dum*rrefoa/(12.E0*dx)
    if(neo_n_radial.eq.1)gradneophi(i)=0.E0
   
    do is = 1, neo_nsp
      do j = 1,nmu
        do k = 1,nvpar
          
          !Second order finite difference radially
          dum = 8*(neof(4,is,i,j,k)-neof(2,is,i,j,k))
          dum = dum + (-neof(5,is,i,j,k)+neof(1,is,i,j,k))
          gradneof(is,i,j,k) = dum*rrefoa/(12.E0*dx)
          if(neo_n_radial.eq.1)gradneof(is,i,j,k)=0.E0

          !In the parallel velocity and parallel directions
          !we use the same scheme as in the code in general  
          !These are naturally more suited to be placed in 
          !linear terms module     

        enddo
      enddo 
    enddo
  enddo

end subroutine cal_gradients

!-------------------------------------------------------------------------
!> Interpolates the data from neo theta grid (which is uniform in theta)
!> onto the GKW theta grid (which is noninform in theta when shaping is used)
!>
!>Interpolation required in the theta direction as the grids
!>do not match up
!-------------------------------------------------------------------------
subroutine neointerp

  use geom,   only : interpquad, pol_angle
  use grid,   only : n_s_grid, nsg_pt

  integer :: is, ix, i, j, k
  real    :: dum(1:nsg_pt), dum2(1:nsg_pt) 
  real    :: coeff_in(1:neo_n_theta+1) 

  do k=1,n_s_grid
    dum2(k) = pol_angle(1,k)
  enddo
  do ix=1,neo_n_radial
    do is=1,neo_n_species
      do i=1,neo_n_energy+1
        do j=1,neo_n_xi+1
          
          !Must give the exact ranges as some arrays have been 
          !extended for parallelisation, or symmetry reasons.
         
          do k = 1,neo_n_theta
            coeff_in(k) = neo_coeffs(ix,is,i,j,k)
          enddo
          coeff_in(neo_n_theta+1)= neo_coeffs(ix,is,i,j,1)        
          
          call interpquad(neo_theta,coeff_in, &
             & neo_n_theta+1,nsg_pt,dum2, &
             & dum)
          
          do k =1,nsg_pt
            neo_coeffs_inter(ix,is,k,i,j)=dum(k)
          enddo
 
          !Ensure periodicity in the s direction (Assuming fourth order)
          neo_coeffs_inter(ix,is,nsg_pt+1,i,j)=neo_coeffs_inter(ix,is,1,i,j)
          neo_coeffs_inter(ix,is,nsg_pt+2,i,j)=neo_coeffs_inter(ix,is,2,i,j)

          neo_coeffs_inter(ix,is,0,i,j)=neo_coeffs_inter(ix,is,nsg_pt,i,j)
          neo_coeffs_inter(ix,is,-1,i,j)=neo_coeffs_inter(ix,is,nsg_pt-1,i,j)
         
        enddo
      enddo
    enddo
  enddo  

  !Then interpolate the electrostatic potential
  do ix=1,neo_n_radial

    call interpquad(neo_theta,neophi_G(ix,:), &
       & neo_n_theta+1,nsg_pt,dum2,dum)

    do k =1,nsg_pt
      neophi_inter(ix,k)=dum(k)
    enddo

    neophi_inter(ix,nsg_pt+1)=neophi_inter(ix,1)
    neophi_inter(ix,nsg_pt+2)=neophi_inter(ix,2)

    neophi_inter(ix,0)=neophi_inter(ix,nsg_pt)
    neophi_inter(ix,-1)=neophi_inter(ix,nsg_pt-1) 
  enddo
  
endsubroutine neointerp

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine neoequil_init
 
  use mpiinterface, only : root_processor

  !Allocates and reads data from NEO output
  call read_neo_gridsizes
  ! using the values read so far, allocate buffers so that the rest can be read, too.
  call neoequil_allocate  
  ! read coefficients for fields from NEO output
  call read_neo_f
  call read_neo_phi
  !Interpolation required in the theta direction as the grids
  !do not match up
  call neointerp  
  !Reconstructs the distribution function and transforms into
  !GKW velocity space
  call reconstruct_f
  !Normalises all distribution fuction etc
  call neof_renormalise

  !Radial gradient is calculated here
  call cal_gradients 
   
  !Deallocate the arrays that are used in reading and transforming.
  call neoequil_deallocate_tmps

  if(root_processor)write(*,*) 'NEO data successfully read and transformed'

end subroutine neoequil_init

end module neoequil
