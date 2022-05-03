!-----------------------------------------------------------------------------
!> Calculates the collision operator coefficients
!> Adds the linear terms of the collision operator
!-----------------------------------------------------------------------------
module collisionop

  use global,       only : lenswitch

  implicit none

  private

  public :: collision_operator_setup,mom_conservation,ene_conservation
  public :: collisionop_read_nml,collisionop_write_nml,collisionop_bcast_nml
  public :: collisionop_check_params,coll_mom_change_int
  public :: coll_mom_change_int_numu,conservation
  public :: cons_type

  interface collisionop_write_nml
    module procedure collisionop_read_nml
  end interface

  !> The reference major radius (only used in the collision operator)
  real, save :: rref

  !> The reference temperature in units of keV (only used for the
  !> collision operator
  real, save :: tref

  !> The reference density in units 10^19 m^-3 (only used for the
  !> collision operator
  real, save :: nref

  !> Scale the electron-ion scattering term to match the effective
  !> ion charge specified zeff specified in the collisions namelist
  !> If impurity species are included,  the input zeff is decreased
  !> accordingly and only electron/main ion collisions are modified
  !> This is distinct from zeff_sp of the input species
  real, save :: zeff

  !Switches for the various terms within the collision operator
  !>Pitch angle scattering
  logical, save :: pitch_angle = .false.
  !> Energy scattering
  logical, save :: en_scatter = .false.
  !> Friction
  logical, save :: friction_coll = .false.
  !> Switch whether or not we desire momentum to be conserved in the collision
  !> operator.  Significant slow down due to integrations.
  logical, save :: mom_conservation
  logical, save :: ene_conservation
  logical, save :: conservation

  !> Switch for the type of ad-hoc field particle collision operator:
  !> - type 'Xu' described in Xu and Rosenbluth Phys Fluids B, 627 (1991)
  !> - type 'Lin' described in Z. Lin et al Physics of Plasmas 2, 2975 (1995)
  character (len = lenswitch), save :: cons_type
  
  !> Switch that reduces the collision operator to a Lorentz operator. i.e.
  !> Just pitch-angle scattering on the electrons by the background ions.
  logical, save :: lorentz = .false.

  !> Values of various moments over the Maxwellian
  !Zeroth order moment
  real, save :: partnum
  !Parallel velocity moment
  real, save :: parmom
  !v^2 moment
  real, save ::ene
  !v^4 moment
  real, save :: enesqr
  !Velocity space integral of: Fmax*(erf(x)-2x*erfp(x))/x and 
  !Fmax*(erf(x)-x*erfp(x))/x^3 used for cons_type='Lin'
  real, save :: erfinte
  real, save :: erfintm

  !> Mass conservation requires a slightly unphysical boundary condition
  !> the diffusion coefficients on the boundaries are set to zero for zero
  !> flux across the boundary.
  !> If false the boundaries are open and mass can flow out of the domain
  !> and therefore mass is no longer conserved.
  !> Mass conservation only works to machine precision when regular mu spacing
  !> is used
  logical, save :: mass_conserve

  !> Parameter from the collision namelist, set to true of user defined
  !> collision frequency is wanted
  logical, save :: freq_override

  !> User defined collision frequecy, which in our units should be fed in for
  !> a Z=1 species at the reference thermal velocity
  real, save :: coll_freq
   
  
  !> Logical switch for momentum conservation.  If true only self collisions
  !> conserve momentum, otherwise all collisions conserve.
  logical, save :: selfcollcon

  !> The factor Gamma^a/b of the collision operator
  real, save, allocatable :: gammab(:,:,:,:)

  !> debug file handle
  integer, save :: i_col
  
contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This subroutine reads (or writes) the reference value of the mass the
!> density and the temperature. These are needed only to calculate the
!> collision frequency.
!----------------------------------------------------------------------------
subroutine collisionop_read_nml(lun,io_stat,lwrite)
  use io, only : write_run_parameter
  use control, only : lcollisions
  use mpiinterface, only : root_processor

  integer, intent(in)           :: lun
  integer, intent(out)          :: io_stat
  logical, optional, intent(in) :: lwrite

  namelist /collisions/ nref, tref, rref , pitch_angle, en_scatter,          &
      & friction_coll, mom_conservation, ene_conservation, mass_conserve,    &
      & freq_override,  coll_freq, selfcollcon, zeff, lorentz,               &
      & cons_type

  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then

      ! Default values..overidden by namelist read
      rref = 1.
      tref = 1.
      nref = 1.
      zeff = 1.
      pitch_angle = lcollisions
      en_scatter =  lcollisions
      friction_coll = lcollisions
      mom_conservation = .false.
      ene_conservation = .false.
      cons_type = 'Xu'
      mass_conserve = .false.
      freq_override = .false.
      coll_freq =0.
      selfcollcon = .true.
      lorentz = .false.

      read(lun,NML=collisions,IOSTAT=io_stat)
    else
      ! do nothing
    end if
  else
    if(root_processor) write(lun,NML=collisions)

    call write_run_parameter('collisions', 'nref', nref)
    call write_run_parameter('collisions', 'tref', tref)
    call write_run_parameter('collisions', 'rref', rref)
    call write_run_parameter('collisions', 'pitch_angle', pitch_angle)
    call write_run_parameter('collisions', 'en_scatter', en_scatter)
    call write_run_parameter('collisions', 'friction_coll', friction_coll)
    call write_run_parameter('collisions', 'mom_conservation', mom_conservation)
    call write_run_parameter('collisions', 'ene_conservation', ene_conservation)
    call write_run_parameter('collisions', 'cons_type', cons_type)
    call write_run_parameter('collisions', 'mass_conserve', mass_conserve)
    call write_run_parameter('collisions', 'freq_override', freq_override)
    call write_run_parameter('collisions', 'coll_freq', coll_freq)
    call write_run_parameter('collisions', 'selfcollcon', selfcollcon)
    call write_run_parameter('collisions', 'zeff', zeff)
    call write_run_parameter('collisions', 'lorentz', lorentz)

  end if

end subroutine collisionop_read_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> bcast the collisionop namelist params
!----------------------------------------------------------------------------
subroutine collisionop_bcast_nml

  use mpiinterface, only : mpibcast
  
  call mpibcast(rref,                       1)
  call mpibcast(tref,                       1)
  call mpibcast(nref,                       1)
  call mpibcast(zeff,                       1)
  call mpibcast(pitch_angle,                1)
  call mpibcast(en_scatter,                 1)
  call mpibcast(friction_coll,              1)
  call mpibcast(mom_conservation,           1)
  call mpibcast(ene_conservation,           1)
  call mpibcast(cons_type,          lenswitch)
  call mpibcast(mass_conserve,              1)
  call mpibcast(freq_override,              1)
  call mpibcast(coll_freq,                  1)
  call mpibcast(selfcollcon,                1)
  call mpibcast(lorentz,                    1)

end subroutine collisionop_bcast_nml

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> put any checks that can be done before memory allocation in here
!----------------------------------------------------------------------------
subroutine collisionop_check_params

  use control,    only : lcollisions, flux_tube
  use components, only : zeff_sp
  use general,    only : gkw_exit, gkw_warn 
  use global,     only : r_tiny
  
  if (lorentz) then
    pitch_angle=.true.
    en_scatter=.false.
    friction_coll=.false.
  endif

  if ((en_scatter .or. friction_coll .or. pitch_angle)) then
    if (.not. lcollisions) then    
      call gkw_exit('Collisions: set collisions=T also in CONTROL')
    end if
  else  
    if (lcollisions) then
      call gkw_warn('Collisions: set collisions=F also in CONTROL')  
    end if
  end if  

  if (en_scatter .and. .not. friction_coll) then
    call gkw_exit('Collisions: Energy scatter requires friction term')
  end if
  
  if (friction_coll .and. .not. en_scatter) then
    call gkw_exit('Collisions: Friction term requires energy scatter')
  end if 

  if (freq_override) then
    !root_and_verbose not yet (re)set
    call gkw_warn('User defined collision frequencies selected - '//     &
        &         'rref input ignored')
    call gkw_warn('tref, nref still used in Couloumb logarithm')
    !Change values to be written back to input.out
    rref = 0.
  end if

  if (nref < r_tiny) then
    call gkw_exit('User defined nref has to be strictly positive')
  end if

  if (tref < r_tiny) then
    call gkw_exit('User defined tref has to be strictly positive')
  end if

  if (ene_conservation .and. .not.mass_conserve) then
    call gkw_warn('Collisions: Energy conservation is only exact when &
                 & mass_conserve velocity space boundaries are also used')
  end if
 
  if (mom_conservation.or.ene_conservation) then
    conservation=.true.
  else
    conservation=.false.
  end if
  
  if (zeff < 1.0) then
    call gkw_exit('Collisions: zeff<1 not allowed')
  end if

  if (.not. flux_tube) then
    call gkw_warn('Collisions: Zeff scalling not implemented for global runs, &
                 & changed to zeff=1.0 and zeff_sp=1.0')
    zeff = 1.0
    zeff_sp = 1.0
  end if

  if (zeff < zeff_sp) then
    call gkw_warn('Collisions: Species Zeff higher than collisions Zeff &
                 & - no enhancement applied to the electron/main ion pitch angle scattering')
    zeff = zeff_sp
  end if
  ! Instead of the above, one could scale down the scattering term, but this may cause confusion
  ! Just kept here in case somebody is interested...
  ! if (zeff < zeff_sp) then
  !  call gkw_warn('Collisions: Species Zeff higher than collisions Zeff &
  !               & - the electron/main ion pitch angle scattering will be scaled DOWN')
  !end if

end subroutine collisionop_check_params

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Collisions initialisation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine collision_operator_setup
  use grid,       only : number_of_species,ns,nx
  use control,    only : lcollisions,vp_trap,uniform_mu_grid
  use general,    only : gkw_abort
  use components, only : nsps
  use global,     only : root_and_verbose
  use mpiinterface, only : root_processor
  use io,         only : get_free_file_unit, output_enabled

  !Error parameter
  integer ierr
  ierr = 0
  
  if (.not. lcollisions) return
  
  !All the collision operator setup is performed here
  if (.not. allocated(gammab)) then
    allocate(gammab(nx,number_of_species,nsps,ns), stat = ierr)
    if (ierr /= 0) then
      call gkw_abort('Could not allocate gammab in collisions')
    end if
  end if

  if (root_processor .and. output_enabled) then
    ! empty the file.
    call get_free_file_unit(i_col)
    open(UNIT=i_col, FILE='Coll_params.dat', &
       & FORM='formatted', STATUS='replace', POSITION='rewind')
    close(i_col)
  end if

  ! call the routine that initializes the factor for the collision operator
  call collision_init
  
  if(vp_trap.eq.0)then
    
    if(uniform_mu_grid)then    
      call collision_differential_terms 
      if(conservation)then
        call cons_momentum
      end if 
    else
      call collision_differential_numu
      if(conservation)then
        call cons_momentum
      end if
    endif

  else if(vp_trap.eq.1)then
    call gkw_abort('There is no collision operator for nonuniform velocity'//&
        &          ' grids')
    !call collision_op_nonuni
  else
    call gkw_abort('Error in Collision operator call')
  end if

  if(root_processor) close(i_col)
  
  ! APS: this can be deallocated in any case now, right?
  if(.not. conservation)then
    if (allocated(gammab)) deallocate(gammab)
  end if

  if (root_and_verbose) then
    write(*,*)'Collision operator initialisation complete'
    write(*,*)
  end if

end subroutine collision_operator_setup


!--------------------------------------------------------------------
!
!> This subroutine caclulates Gamma(a,b) necessary for the collision 
!> operator. Gamma is defined as 
!>
!>  gammab = R_ref nb q_a^2 q_b^2 lambda^{a/b} /
!>           (4 pi epsilon_0^2 m_a^2 vth_a^4) 
!> 
!> where R_ref is the reference major radius [in meters] 
!>       nb is the density of the species b
!>       q_a is the particle charge of species a (scattered)
!>       q_b is the charge of the particle of species b (scattering)
!>       lambda is the Coulomb logarithm
!>       m_a is the particle mass (species a)  
!>       vth_a is the thermal velocity (species a) 
!
!--------------------------------------------------------------------
subroutine collision_init 

  use grid,         only : number_of_species, ns, nx
  use components,   only : de_G,mas_G,tmp_G,signz_G,nsps, zeff_sp
  use rotation,     only : cfen_G, rotation_parallelized
  use io,           only : open_real_lu, close_lu, output_enabled
  use io,           only : append_chunk, xy_fmt, ascii_fmt, attach_metadata
  use general,      only : gkw_abort, gkw_warn
  use mpiinterface, only : root_processor
  use global,       only : root_and_verbose
  
  ! integers for the loops 
  integer :: i, j, k, ix, ierr, nxmax

  ! index for the main ion species
  integer :: i_msp

  ! pitch angle scattering enhancement for electron/main ion collisions
  real :: zeff_pa

  ! the electron density 
  real :: iigamma
  real, allocatable :: electron_density(:)

  !The reference gamma calculated from the input collision frequency
  real :: ref_freq

  integer :: lun

  if(.not.rotation_parallelized) then
    call gkw_abort('collision init must be called after parallelize rotation')
  end if

  allocate(electron_density(nx), stat = ierr) 
  if (ierr /= 0) call gkw_abort('Can not allocate electron_density')

  ! determine the electron density in units of 10^19 
  ! (used only for the coulomb logarithm)
  ! find the electron species
  electron_density(:) = 0.
  do j = 1, nsps
    if (signz_G(j) < 0) then 
      do ix = 1, nx 
        electron_density(ix) = electron_density(ix) + de_G(ix,j) * nref
      end do   
    endif
  end do

  ! calculate the Coulomb logarithm. Note the calculation is stored 
  ! in gammab, which is later transformed in the pre-factor gamma(a,b)
  ! of the collision operator 

  ! i is the scattered species
  ! j is the scattering species
  do ix = 1, nx 
    do i = 1, number_of_species
      do j = 1, nsps
        ! check for electrons 
        if (signz_G(i) < 0.) then
          if (signz_G(j) < 0.) then 
            ! electron electron collisions
            !The factor of 0.1 in the expression is to account for densities
            !being given in 10^19 while in the literature they are 10^20 
            if(lorentz)then            
              gammab(ix,i,j,:) = 0.E0
            else
              gammab(ix,i,j,:) = 14.9E0 - 0.5E0*log(0.1E0*electron_density(ix)) +  &
                & log(tmp_G(ix,i)*tref)
            endif
          else 
            ! electron ion collisions (from NRL p.34)
            if((tmp_G(ix,i)*tref).lt.(0.01*signz_G(j)**2))then
              gammab(ix,i,j,:) = 17.2E0 - 0.5E0*log(0.1E0*signz_G(j)**2*electron_density(ix)) +  &
                 & 1.5*log(tmp_G(ix,i)*tref)
            else
              gammab(ix,i,j,:) = 14.8E0 - 0.5E0*log(0.1E0*electron_density(ix)) +  &
                 & log(tmp_G(ix,i)*tref)
            endif 
          endif
        else  
          ! ion electron collisions  (from NRL p.34)
          if (signz_G(j) < 0.) then
            if((tmp_G(ix,j)*tref).lt.(0.01*signz_G(i)**2))then
              if(lorentz)then            
                gammab(ix,i,j,:) = 0.E0
              else
                gammab(ix,i,j,:) = 17.2E0 - 0.5E0*log(0.1E0*signz_G(i)**2*electron_density(ix)) +  &
                  & 1.5*log(tmp_G(ix,j)*tref)
              endif
            else
              if(lorentz)then            
                gammab(ix,i,j,:) = 0.E0
              else
                gammab(ix,i,j,:) = 14.8E0 - 0.5E0*log(0.1E0*electron_density(ix)) +  &
                  & log(tmp_G(ix,j)*tref)
              endif 
            end if
          else
            ! ion -ion collisions 
            !Old expression
            !gammab(i,j,:) = 17.3E0 - 0.5E0*log(0.1E0*signz_G(i)*signz_G(j)*electron_density) +  &
            !    & 1.5*log(0.5E0*(tmp_G(i)+tmp_G(j))*tref)
            !The expression below is the full NRL formula for ion-ion with arbitrary mass and
            !charge.
            if (de_G(ix,i) > 0. .or. de_G(ix,j) > 0.) then
              if(lorentz)then            
                gammab(ix,i,j,:) = 0.E0
              else
                gammab(ix,i,j,:) = 17.3E0 - 1.E0*log(signz_G(i)*signz_G(j)*(mas_G(i)+mas_G(j))) +  & 
                  & 1.E0*log((mas_G(i)*tmp_G(ix,j)+mas_G(j)*tmp_G(ix,i))*tref) - &
                  & 0.5E0*log(0.1E0*nref/tref) - &
                  & 0.5E0*log((de_G(ix,i)*signz_G(i)**2/tmp_G(ix,i))+(de_G(ix,j)*signz_G(j)**2/tmp_G(ix,j)))
              endif            
            else !Collisions of two zero density species
              gammab(ix,i,j,:) = 0.0
            endif
          endif
        endif
      end do
    end do
  end do 

  if (root_and_verbose) then
    write(*,*)'!Coulomb logarithms!'
    write(*,*)'!First index scattered species'
    write(*,*)'!Second index is the scattering species!'
  end if
  if (root_and_verbose) then
    !do ix = 1, nx 
    !write(*,*)'Grid point ',ix 
    ix = 1
    do i = 1, number_of_species
      do j = 1, nsps
        write(*,*) i,j,gammab(ix,i,j,1)
      end do
    end do
    !end do
  end if

  !Write collision parameters to a file for future reference
  nxmax=1;
  !If need to output Coll_parmas for all x, consider
  !a simple solution to parallel_x output issue using proc_subset;
  !See how it is used in collision_differential_terms
  !nxmax = n_x_grid; if (flux_tube) nxmax = 1
  if (root_processor .and. output_enabled) then
    ! FIXME What is the correct STATUS setting?
    open(UNIT=i_col, FILE='Coll_params.dat', &
         & FORM='formatted', POSITION='append')

    write(i_col,*)'Coulomb Logarithms: scattered, scattering, value'
    call open_real_lu('CoulombLog', 'collisionop', (/ 3/), ascii_fmt, lun)
    call attach_metadata(lun, 'description', &
       & 'Coulomb Logarithms: scattered, scattering, value', ascii_fmt)
    do i = 1, number_of_species 
      do ix = 1, nxmax; do j = 1, nsps 
        write(i_col,*)i,j,gammab(ix,i,j,1)
        call append_chunk(lun, (/ i*1.0, j*1.0, gammab(ix,i,j,1) /), &
           & xy_fmt, ascii_fmt)
      end do; end do 
    end do
    call close_lu(lun, ascii_fmt)
  end if 

  if(freq_override)then
    ! The input collision frequency is assumed to be the ion-ion
    ! scattering freqency.  Therefore the background density is the ion
    ! density and the scattered species mass is the ion mass.
    if((signz_G(1).lt.1.E0).or.(signz_G(1).gt.1.E0))then
       write(*,*) 'Charge of first species', signz_G(1)
       call gkw_abort('When user defined collision frequencies are chosen, &
         & the first species must be singly positively charged')
    end if

    iigamma =  gammab(1,1,1,1)
    ref_freq = coll_freq
    do k=1,ns
      do i=1,number_of_species
        do j=1,nsps; do ix = 1, nx 
          ! Input freqency (single charged ion-ion frequency) should be
          ! calculated from
          !   nu = 6.5141x10^{-5}*Rref*logLii*nref/Tref^2
          ! and fed in as ref_freq.  This is the ion-ion frequency, then is
          ! rescaled according to the expression below.
          gammab(ix,i,j,k) = signz_G(i)**2*signz_G(j)**2*ref_freq*de_G(ix,j)*      &
              &  exp(-cfen_G(k,j))*(gammab(ix,i,j,k)/iigamma)/(1.E0*tmp_G(ix,i)**2)
        end do; end do 
      end do
    end do
  else !not freq_override
    do k=1,ns
      do i=1,number_of_species
        do j=1,nsps; do ix = 1, nx 
          ! Alternatively if freq_override is false it is calculated as
          ! follows, which is kinda the same.
          gammab(ix,i,j,k) = 6.5141e-5*(rref*nref/(tref**2))*de_G(ix,j)      &
              & *exp(-cfen_G(k,j))*signz_G(i)**2*signz_G(j)**2               &
              & *gammab(ix,i,j,k)/(tmp_G(ix,i)**2)
        end do; end do 
      end do
    end do
  end if !freq_override
 
  ! scale the electron-main ion scattering by additional fictional zeff
  ! (assumes electrons have negative charge)
  ! The scaling factor is calculated as
  ! zeff_pa = 1 + (zeff - zeff_sp)*ne/(ni*Zi)
  ! with i the main ion species index 
  ! Not made for GLOBAL runs yet 
  i_msp    = maxloc(de_G(1,:),1,(signz_G(:)>0))  ! index of main species ion
  zeff_pa = 0.0
  do i=1,number_of_species
   if (signz_G(i) < 0 ) then                 ! electron species
     zeff_pa = zeff_pa + de_G(1,i)*(zeff-zeff_sp)
   end if  
  end do
  zeff_pa = 1.0 + zeff_pa / (de_G(1,i_msp)*signz_G(i_msp)**2)
  if ( zeff_pa < 0) then
   call gkw_warn('Collisions: zeff_pa set to zero in e/i pitch-angle scattering &
                &scaling to avoid negative collision frequency')
  end if

  do i=1,number_of_species; do j=1,nsps
    if (signz_G(i) < 0 .and. j/=i .and. j==i_msp) then
      do ix = 1, nx  
        gammab(ix,i,j,:)=zeff_pa*gammab(ix,i,j,:)
      end do 
    end if
  end do; end do

  if (root_processor) then
     write(*,*)
     write(*,*) 'Electron/main-ion collision frequency &
                &scattering scaled up by Zeff=', zeff_pa
     write(*,*)
  end if

  if (root_and_verbose) then
      write(*,*)'!Normalised collision frequencies!'
      write(*,*)'!First index scattered species'
      write(*,*)'!Second index is the scattering species!'
  end if
  
  if (root_processor) then
    write(i_col,*)'Collision Frequencies: scattered, scattering, Freq (Vth(sp)/Rref)'
    call open_real_lu('CollFreqs', 'collisionop',(/ 3/), ascii_fmt, lun)
    call attach_metadata(lun, 'description_2', &
       & 'Collision Frequencies: scattered, scattering, Freq (Vth(sp)/Rref)', ascii_fmt)
    do i = 1, number_of_species
      do ix = 1, nxmax ; do j = 1, nsps
        write(i_col,*)i,j,gammab(ix,i,j,1)
        call append_chunk(lun, (/ i*1.0, j*1.0, gammab(ix,i,j,1) /), &
           & xy_fmt, ascii_fmt)
      end do; end do
    end do
    call close_lu(lun, ascii_fmt)
    close(i_col)
  end if

  ix = 1
  do i = 1, number_of_species
    do j = 1, nsps  
      if (root_and_verbose) write(*,*)i,j,gammab(ix,i,j,1)
      end do 
  end do  

  ! deallocate help array 
  deallocate(electron_density) 

end subroutine collision_init

!--------------------------------------------------------------------
!> This routine calculates the diffusion coefficient of pitch angle 
!> scattering. Note dthth is not exactly D_theta,theta one usually 
!> finds in the literature. It is defined here as 
! 
!> dthth = sum_b (1/4) gamma^(a/b)_N [ (2-1/vtb**2) erf(vtb) + 
!         erf^prime(vtb) / vtb ] 
!
!--------------------------------------------------------------------
subroutine caldthth(vp,mubn,ix,is,dthth,k)

  use components, only : vthrat, vthrat_G, nsps
  use grid,       only : gsp
  use specfun,    only : erf => sf_erf

  real,  intent(in)   :: vp, mubn 
  integer, intent(in) :: ix,is, k  
  real, intent(out)   :: dthth 

  integer :: i 
  real vn, vtb, dum

  vn = sqrt(vp**2 + 2.E0*sqrt(mubn**2)) 

  dthth = 0. 
  ! N.B. _is_ corresponds to the local species
  do i = 1, nsps
    vtb = vn*vthrat(is) / vthrat_G(i) 
    dum = (2.E0 - 1.E0/(vtb**2))*erf(vtb) + erfp(vtb)/vtb 
    dthth = dthth + gammab(ix,gsp(is),i,k)*dum/(4.E0*vn) 
  end do

end subroutine caldthth


!--------------------------------------------------------------------
!> This routine calculates the diffusion coefficient of energy 
!> scattering. It is defined here as 
! 
!> dvv = sum_b (1/2) gamma^(a/b)_N [ (1/vtb**2) erf(vtb) - 
!         erf^prime(vtb) / vtb ] 
!
!--------------------------------------------------------------------
subroutine caldvv(vp,mubn,ix,is,dvv,k)

  use components, only : vthrat, vthrat_G, nsps
  use grid,       only : gsp
  use specfun,    only : erf => sf_erf

  real,  intent(in)   :: vp, mubn 
  integer, intent(in) :: ix,is,k  
  real, intent(out)   :: dvv 

  integer i 
  real vn,  vtb, dum

  vn = sqrt(vp**2 + 2.E0*mubn) 

  dvv = 0. 
  do i = 1, nsps 
    vtb = vn * vthrat(is) / vthrat_G(i) 
    dum = erf(vtb)/(vtb**2) - erfp(vtb)/vtb 
    dvv = dvv + gammab(ix,gsp(is),i,k)*dum/(2.E0*vn) 
  end do

end subroutine caldvv


!--------------------------------------------------------------------
!> This routine calculates the diffusion coefficient of energy 
!> scattering. It is defined here as 
! 
!> dvv = sum_b (1/2) gamma^(a/b)_N [ (1/vtb**2) erf(vtb) - 
!         erf^prime(vtb) / vtb ] 
!
!--------------------------------------------------------------------
subroutine selfcaldvv(vp,mubn,ix,is,dvv,k)

  use grid,       only : gsp
  use specfun,    only : erf => sf_erf

  real,  intent(in)   :: vp, mubn 
  integer, intent(in) :: ix,is, k  
  real, intent(out)   :: dvv 

  real vn, dum

  vn = sqrt(vp**2 + 2.E0*mubn) 

  dvv = 0.  
  dum = erf(vn)/(vn**2) - erfp(vn)/vn 
  dvv = gammab(ix,gsp(is),gsp(is),k)*dum/(2.E0*vn) 
  
end subroutine selfcaldvv


!--------------------------------------------------------------------
!> This is a copy of the above caldthth function, but instead of sum
!> ming over all species only considers self collisions for the 
!> momentum conserving term
!--------------------------------------------------------------------
subroutine selfcaldthth(vp,mubn,ix,is,dthth,k)

  use grid,       only : gsp
  use specfun,    only : erf => sf_erf

  real,  intent(in) :: vp, mubn 
  integer, intent(in) :: ix,is,k  
  real   dthth 

  real vn, dum 

  vn = sqrt(vp**2 + 2.E0*sqrt(mubn**2)) 
  dthth = 0.   
  dum = (2.E0 - 1.E0/(vn**2))*erf(vn) + erfp(vn)/vn 
  dthth = gammab(ix,gsp(is),gsp(is),k)*dum/(4.E0*vn) 
  
end subroutine selfcaldthth

!!$!--------------------------------------------------------------------
!!$!> This routine calculates the friction coefficient. It is defined 
!!$!> here as 
!!$! 
!!$!> fv = -sum_b (ma/mb) gamma^(a/b)_N [ erf(vtb) - 
!!$!         vtb erf^prime(vtb)  ] 
!!$!
!!$!--------------------------------------------------------------------
subroutine calfv(vp,mubn,ix,is,fv,k)

  use components, only : vthrat, mas, mas_G, vthrat_G, nsps
  use grid,       only : gsp
  use specfun,    only : erf => sf_erf
  
  real,  intent(in) :: vp, mubn 
  integer, intent(in) :: ix,is,k 
  real   fv,mrat

  integer i 
  real vn,  vtb, dum

  !Need to calculate the mass ratio.

  vn = sqrt(vp**2 + 2.E0*mubn) 

  fv = 0. 
  do i = 1, nsps
    mrat = mas(is)/mas_G(i)
    vtb = vn * vthrat(is)/vthrat_G(i)
    dum = erf(vtb) - erfp(vtb)*vtb 
    fv = fv + mrat*gammab(ix,gsp(is),i,k)*dum/(vn**2)    
  end do

end subroutine calfv


!--------------------------------------------------------------------
!> This routine calculates the friction coefficient. It is defined 
!> here as 
!> 
!> fv = -sum_b (ma/mb) gamma^(a/b)_N [ erf(vtb) - 
!>         vtb erf^prime(vtb)  ] 
!>
!--------------------------------------------------------------------
subroutine selfcalfv(vp,mubn,ix,is,fv,k)

  use grid,       only : gsp
  use specfun,    only : erf => sf_erf

  real,  intent(in)   :: vp, mubn 
  integer, intent(in) :: ix, is,k
  real, intent(out)   :: fv

  real vn, dum

  !Need to calculate the mass ratio.

  vn = sqrt(vp**2 + 2.E0*mubn) 

  fv = 0. 
  dum = erf(vn) - erfp(vn)*vn 
  fv = fv + gammab(ix,gsp(is),gsp(is),k)*dum/(vn**2)  
  
end subroutine selfcalfv

!--------------------------------------------------------------------
!> output coefficients for every species
!--------------------------------------------------------------------
subroutine output_coefficients
  use mpiinterface, only : mpibarrier, mpiallreduce_sum_inplace
  use grid, only : lsp, number_of_species, proc_subset
  use io, only : open_real_lu, close_lu, append_chunk, xy_fmt
  use io, only : ascii_fmt, get_free_file_unit
  integer :: is, lun
  real :: Dtemp
  do is = 1, number_of_species
    call mpibarrier

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      ! call get_free_file_unit(i_col)
      ! ! FIXME What is the correct STATUS setting?
      ! open(UNIT=i_col, FILE='Coll_params.dat', &
      !    & FORM='formatted', POSITION='append')

      call caldthth(1.E0,0.E0,1,lsp(is),Dtemp,1)
      ! write(i_col,*)'D_theta_theta'
      ! write(i_col,*) is,Dtemp
    end if
    call open_real_lu('D_theta_theta', 'collisionop', (/ 2/), ascii_fmt, lun)
    call mpiallreduce_sum_inplace(Dtemp,1) ! ~ send scalar to root process
    call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
    call close_lu(lun, ascii_fmt)

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      call caldvv(1.E0,0.E0,1,lsp(is),Dtemp,1)
      ! write(i_col,*)'D_nu_nu'
      ! write(i_col,*) is,Dtemp
    end if
    call open_real_lu('D_nu_nu', 'collisionop', (/ 2/), ascii_fmt, lun)
    call mpiallreduce_sum_inplace(Dtemp,1) ! ~ send scalar to root process
    call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
    call close_lu(lun, ascii_fmt)

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      call calfv(1.E0,0.E0,1,lsp(is),Dtemp,1)
      ! write(i_col,*)'F_nu'
      ! write(i_col,*) is,Dtemp
      ! close(i_col)
    end if
    call open_real_lu('F_nu', 'collisionop', (/ 2/), ascii_fmt, lun)
    call mpiallreduce_sum_inplace(Dtemp,1) ! ~ send scalar to root process
    call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
    call close_lu(lun, ascii_fmt)

  end do
end subroutine output_coefficients


!----------------------------------------------------------------
!> Collision operator linear terms
!> This version of the operator assumes that a uniform mu grid 
!> has been used. 
!----------------------------------------------------------------
subroutine collision_differential_terms

  use control,        only : neoclassics
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp, gvpar
  use grid,           only : gmu, n_vpar_grid, n_mu_grid
  use grid,           only : lsp, proc_subset
  use index_function, only : indx
  use geom,           only : bn
  use components,     only : vthrat
  use matdat,         only : set_indx, add_element
  use velocitygrid,   only : mugr, vpgr, dvp, dmu, vgridboundary
  use dist,           only : ifdis
  use structures,     only : matrix_element

  integer :: imod, ix, i, j, k, is, ierr, jjh
  logical :: ingrida
  character(len=64) :: term='Collisions differential'
  type (matrix_element) :: elem

  !Values of diffusion coefficient in the 4 half grid points of interest
  !and the point of consideration
  !Temporary storage of the diffusion coefficients and the prefactors
  real :: Dtemp
  !Half point velocity values
  real :: vdum, mudum
  real :: fac,fad,faf,vzero,del
  complex :: mat_elem

  ! Identifier of the term 
  elem%term = term
  elem%ideriv = 2

  ! Type of iih and jjh 
  elem%itype = ifdis 
  elem%itloc = ifdis

  call output_coefficients
  
  do is = 1,nsp   
    do imod = 1, nmod
      do ix = 1,nx
        do i = 1,ns
          do j = 1,nmu
            do k = 1,nvpar

              call set_indx(elem,imod,ix,i,j,k,is)

              vzero = sqrt(vpgr(i,j,k,is)**2 + 2.E0*mugr(j)*bn(ix,i))
              
              !First is the dvpardvpar term            
              vdum = vpgr(i,j,k,is)+0.5E0*dvp
              mudum = mugr(j)
              if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = 2.E0*mudum*Dtemp*bn(ix,i)       & 
                     & /(vdum**2 + 2.E0*mudum*bn(ix,i))
                  fac = fac/(dvp*dvp)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fad = vdum**2*Dtemp/sqrt(vdum**2+2*mudum*bn(ix,i)) 
                  fad = fad/(dvp*dvp)                  
                else
                  fad = 0.E0
                endif
                
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  faf = vdum*Dtemp/(dvp)
                else
                  faf = 0.E0
                endif

              endif
 
              call ccdelta(fac+fad,faf,dvp,del)
              
              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad+(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad+del*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if
     
              vdum = vpgr(i,j,k,is)-0.5*dvp
              mudum = mugr(j)     
              if((gvpar(k).eq.1).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = 2.E0*mudum*Dtemp*bn(ix,i)       & 
                     & /(vdum**2 + 2*mudum*bn(ix,i))
                  fac = fac/(dvp*dvp)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fad = vdum**2*Dtemp/sqrt(vdum**2+2*mudum*bn(ix,i)) 
                  fad = fad/(dvp*dvp)     
                else
                  fad = 0.E0
                endif
                
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  faf = vdum*Dtemp/(dvp)     
                else
                  faf = 0.E0
                endif

              end if
          
              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
               if(neoclassics) then
                 call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              !Second is the dmudmu term    
           
              vdum = vpgr(i,j,k,is)
              mudum= mugr(j)+0.5E0*dmu
              if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)    
                  fac = 2.E0*mudum*Dtemp*vdum**2     & 
                      & /((vdum**2 + 2.E0*mudum*bn(ix,i))*bn(ix,i))
                  fac = fac/(dmu*dmu)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fad = 4.E0*mudum**2*Dtemp/sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                  fad = fad/(dmu*dmu)
                else
                  fad = 0.E0
                end if
                  
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  faf = 2.E0*mudum*Dtemp/(dmu)
                else
                  faf = 0.E0
                endif

              endif
              
              call ccdelta(fac+fad,faf,dmu,del)
      
              elem%jloc=j+1
              elem%kloc=k               
              call vgridboundary(1,elem,ingrida) 
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad+ (1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad+del*faf) 
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              vdum = vpgr(i,j,k,is)
              mudum = mugr(j)-0.5E0*dmu
              
              !No need for imposition of mass-conservation on the lower
              !mu boundary as the flux is naturally zero at this point.
              if(pitch_angle)then
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                fac = 2.E0*mudum*Dtemp*vdum**2     &
                   & /((vdum**2 + 2.E0*mudum*bn(ix,i))*bn(ix,i))
                fac = fac/(dmu*dmu)
              else
                fac = 0.E0
              endif  
                
              if(en_scatter)then
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                fad = 4.E0*mudum**2*Dtemp/sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                fad = fad/(dmu*dmu)
              else
                fad = 0.E0
              endif
                
              if(friction_coll)then
                call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                faf = 2.E0*mudum*Dtemp/(dmu)
              else
                faf = 0.E0
              endif

              call ccdelta(fac+fad,faf,dmu,del)
     
              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then               
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k   
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if
              
              !Cross terms (there are no cross terms in friction)
              !Must interpolate (using Bilinear) the half/half grid
              !points
              
              !Firstly dmudvpar
              vdum = vpgr(i,j,k,is)
              mudum = mugr(j)+0.5*dmu
              if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -2*mudum*Dtemp*vdum     & 
                    & /(vdum**2 + 2.E0*mudum*bn(ix,i))
                 else
                   fac = 0.E0
                 endif
                 
                 if(en_scatter)then
                   call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)              
                   fad = 2.E0*mudum*vdum*Dtemp &
                      &  /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                 else
                   fad = 0.E0
                 end if

              end if
              
              elem%jloc=j+1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)                  
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)                     
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if


              vdum = vpgr(i,j,k,is)
              mudum = mugr(j)-0.5E0*dmu       
              if((gmu(j).eq.1).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fac = -2*mudum*vdum*Dtemp     & 
                    & /(vdum**2 + 2.E0*mudum*bn(ix,i))
                else
                  fac = 0.E0
                endif
                 
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)              
                  fad = 2.E0*mudum*vdum*Dtemp &
                      &  /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                else
                  fad = 0.E0
                end if

              end if
            
              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)  
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh =indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              !Then dvpardmu
              vdum = vpgr(i,j,k,is)+0.5E0*dvp
              mudum = mugr(j)
              if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -2.E0*mudum*Dtemp*vdum     & 
                     & /(vdum**2 + 2.E0*mudum*bn(ix,i))
                 else
                   fac = 0.E0
                 endif
                 
                 if(en_scatter)then
                   call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                   fad = 2.E0*mudum*vdum*Dtemp &
                      & /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                 else
                   fad = 0.E0
                 end if
               end if
    
              elem%jloc=j+1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)     
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)    
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              vdum = vpgr(i,j,k,is)-0.5*dvp
              mudum = mugr(j)
              if((gvpar(k).eq.1).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -2*mudum*Dtemp*vdum     & 
                     & /(vdum**2 + 2*mudum*bn(ix,i))
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fad = 2.E0*mudum*vdum*Dtemp &
                      & /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
                else
                  fad = 0.E0
                end if
              end if
              
              elem%jloc=j+1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)       
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)   
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
                end if
              end if

            end do
          end do
        end do

      end do
    end do
  end do

end subroutine collision_differential_terms

!----------------------------------------------------------------
!> Collision operator linear terms
!> This version of the operator assumes that a uniform mu grid 
!> has been used. 
!---------------------------------------------------------------
subroutine collision_differential_numu

  use control,        only : neoclassics
  use grid,           only : nmod, nx, ns, nmu, nvpar, nsp
  use grid,           only : n_vpar_grid, gvpar, n_mu_grid, gmu
  use grid,           only : lsp, number_of_species, proc_subset
  use index_function, only : indx
  use geom,           only : bn
  use components,     only : vthrat
  use matdat,         only : add_element, set_indx
  use velocitygrid,   only : mugr, vpgr, dvp, dvperp, vgridboundary
  use dist,           only : ifdis
  use structures,     only : matrix_element
  use mpiinterface,   only : mpibarrier, mpireduce_sum_inplace, root_processor
  use io,             only : open_real_lu, close_lu, append_chunk, xy_fmt
  use io,             only : ascii_fmt, get_free_file_unit, output_enabled

  integer :: imod, ix, i, j, k, is, ierr, jjh
  logical :: ingrida
  character(len=64) :: term='Collisions differential numu'
  type (matrix_element) :: elem

  !Values of diffusion coefficent in the 4 half grid points of interest
  !and the point of consideration
  !Temporary storage of the diffusion coefficients and the prefactors
  real :: Dtemp
  !Half point velocity values
  real :: vdum, mudum, vpdum
  real :: fac,fad,faf,vzero,del, dvrp
  complex :: mat_elem
  !For both uniform and non-uniform grids
  integer :: lun

  ! Identifier of the term 
  elem%term = term
  elem%ideriv =2 

  ! Type of the terms 
  elem%itype = ifdis 
  elem%itloc = ifdis

  ! output coefficients for every species
  do is = 1, number_of_species
    call mpibarrier

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      call caldthth(1.E0,0.E0,1,lsp(is),Dtemp,1)
      if(output_enabled) then
        ! FIXME What is the correct STATUS setting?
        call get_free_file_unit(i_col)
        open(UNIT=i_col, FILE='Coll_params.dat', &
           & FORM='formatted', POSITION='append')
        write(i_col,*)'D_theta_theta'
        write(i_col,*) is,Dtemp
      end if
    end if
    call mpireduce_sum_inplace(Dtemp) ! ~ send scalar to root process
    if(root_processor) then
      call open_real_lu('D_theta_theta', 'collisionop', (/ 2/), ascii_fmt, lun)
      call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
      call close_lu(lun, ascii_fmt)
    end if

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      call caldvv(1.E0,0.E0,1,lsp(is),Dtemp,1)
      if(output_enabled) then
        write(i_col,*)'D_nu_nu'
        write(i_col,*) is,Dtemp
      end if
    end if
    call mpireduce_sum_inplace(Dtemp) ! ~ send scalar to root process
    if(root_processor) then
      call open_real_lu('D_nu_nu', 'collisionop', (/ 2/), ascii_fmt, lun)
      call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
      call close_lu(lun, ascii_fmt)
    end if

    Dtemp =  0
    if (proc_subset(1,1,1,1,is)) then
      call calfv(1.E0,0.E0,1,lsp(is),Dtemp,1)
      if(output_enabled) then
        write(i_col,*)'F_nu'
        write(i_col,*) is,Dtemp
        close(i_col)
      end if
    end if
    call mpireduce_sum_inplace(Dtemp) ! ~ send scalar to root process
    if(root_processor) then
      call open_real_lu('F_nu', 'collisionop', (/ 2/), ascii_fmt, lun)
      call append_chunk(lun, (/ is*1.0, Dtemp /), xy_fmt, ascii_fmt)
      call close_lu(lun, ascii_fmt)
    end if
  end do

  do is = 1,nsp   
    do imod = 1, nmod
      do ix = 1,nx
        do i = 1,ns
          do j = 1,nmu
            do k = 1,nvpar

              call set_indx(elem,imod,ix,i,j,k,is)
              ! mark the term as a second deriviative (for timestep est)
              elem%ideriv = 2

              vzero = sqrt(2.E0*mugr(j)*bn(ix,i))
              dvrp = sqrt(bn(ix,i))*dvperp
              
              !First is the dvpardvpar term            
              vdum = vpgr(i,j,k,is)+0.5E0*dvp
              mudum = mugr(j)
              vpdum = sqrt(2.E0*bn(ix,i)*mudum)
              if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = vpdum**2*Dtemp       & 
                     & /(vdum**2 + vpdum**2)
                  fac = fac/(dvp*dvp)
                else
                  fac=0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fad = vdum**2*Dtemp/(vdum**2+vpdum**2) 
                  fad = fad/(dvp*dvp)                  
                else
                  fad=0.E0
                endif
                
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  faf = vdum*Dtemp/sqrt(vdum**2 + vpdum**2)
                  faf = faf/dvp
                else
                  faf=0.E0
                endif

              endif
 
              call ccdelta(fac+fad,faf,dvp,del)
              
              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad+(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad+del*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if
     
              vdum = vpgr(i,j,k,is)-0.5*dvp
              mudum = mugr(j)    
              vpdum = sqrt(2.E0*bn(ix,i)*mudum) 
              if((gvpar(k).eq.1).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = vpdum**2*Dtemp       & 
                     & /(vdum**2 + vpdum**2)
                  fac = fac/(dvp*dvp)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fad = vdum**2*Dtemp/(vdum**2+vpdum**2) 
                  fad = fad/(dvp*dvp)     
                else
                  fad = 0.E0
                endif
                
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  faf = vdum*Dtemp/sqrt(vdum**2 + vpdum**2)
                  faf = faf/dvp     
                else
                  faf = 0.E0
                endif

              end if
          
              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              !Second is the dvperpdvperp term    
           
              vdum = vpgr(i,j,k,is)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) + 0.5E0*dvrp
              mudum = vpdum**2/(2*bn(ix,i))
              
              if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
                fac = 0.E0
                fad = 0.E0
                faf = 0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)    
                  fac = vpdum*Dtemp*vdum**2     & 
                      & /(vdum**2 + vpdum**2)
                  fac = fac/(vzero*dvrp*dvrp)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fad = vpdum**3*Dtemp/(vdum**2 + vpdum**2)
                  fad = fad/(vzero*dvrp*dvrp)
                else
                  fad = 0.E0
                end if
                  
                if(friction_coll)then
                  call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  faf = vpdum**2*Dtemp/sqrt(vdum**2 + vpdum**2)
                  faf = faf/(dvrp*vzero)
                else
                  faf = 0.E0
                endif

              endif
              
              call ccdelta(fac+fad,faf,dvrp,del)
      
              elem%jloc=j+1
              elem%kloc=k               
              call vgridboundary(1,elem,ingrida) 
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad+ (1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad+del*faf) 
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if
      
              vdum = vpgr(i,j,k,is)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) - 0.5E0*dvrp
              mudum = vpdum**2/(2*bn(ix,i))
              
              !No need for imposition of mass-conservation on the lower
              !mu boundary as the flux is naturally zero at this point.
              if(pitch_angle)then
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                fac = vpdum*Dtemp*vdum**2     &
                   & /(vdum**2 + vpdum**2)
                fac = fac/(vzero*dvrp*dvrp)
              else
                fac = 0.E0
              endif  
                
              if(en_scatter)then
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                fad = vpdum**3*Dtemp/(vdum**2 + vpdum**2)
                fad = fad/(vzero*dvrp*dvrp)
              else
                fad = 0.E0
              endif
                
              if(friction_coll)then
                call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                faf = vpdum**2*Dtemp/sqrt(vdum**2 + vpdum**2)
                faf = faf/(dvrp*vzero)
              else
                faf = 0.E0
              endif

              call ccdelta(fac+fad,faf,dvrp,del)

              elem%jloc=j
              elem%kloc=k
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then               
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                   call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k   
              call vgridboundary(1,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              !Cross terms (there are no cross terms in friction)
              !Must interpolate (using Bilinear) the half/half grid
              !points
              
              !Firstly dvperpdvpar
              vdum = vpgr(i,j,k,is)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) + 0.5E0*dvrp
              mudum = vpdum**2/(2.E0*bn(ix,i))
              if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -vpdum**2*Dtemp*vdum     & 
                    & /(vzero*(vdum**2 + vpdum**2))
                 else
                   fac = 0.E0
                 endif
                 
                 if(en_scatter)then
                   call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)              
                   fad = vpdum**2*vdum*Dtemp &
                      &  /(vzero*(vdum**2 + vpdum**2))
                 else
                   fad = 0.E0
                 end if

              end if
              
              elem%jloc=j+1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)                  
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)                     
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if


              vdum = vpgr(i,j,k,is)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) - 0.5E0*dvrp
              mudum = vpdum**2/(2*bn(ix,i))    
              if((gmu(j).eq.1).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fac = -vpdum**2*vdum*Dtemp     & 
                    & /(vzero*(vdum**2 + vpdum**2))
                else
                  fac = 0.E0
                endif
                 
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)              
                  fad = vpdum**2*vdum*Dtemp &
                      &  /(vzero*(vdum**2 + vpdum**2))
                else
                  fad = 0.E0
                end if

              end if
            
              elem%jloc=j
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)  
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh =indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                        call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              !Then dvpardvperp
              vdum = vpgr(i,j,k,is)+0.5E0*dvp
              mudum = mugr(j)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j))
              if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -vpdum*vdum*Dtemp     & 
                     & /(vdum**2 + vpdum**2)
                 else
                   fac = 0.E0
                 endif
                 
                 if(en_scatter)then
                   call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                   fad = vpdum*vdum*Dtemp &
                      & /(vdum**2 + vpdum**2)
                 else
                   fad = 0.E0
                 end if
               end if
    
              elem%jloc=j+1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k+1
              call vgridboundary(2,elem,ingrida)     
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)    
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              vdum = vpgr(i,j,k,is)-0.5*dvp
              mudum = mugr(j)
              vpdum = sqrt(2.E0*bn(ix,i)*mugr(j))
              if((gvpar(k).eq.1).and.mass_conserve)then
                fac=0.E0
                fad=0.E0
              else
                if(pitch_angle)then
                  call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
                  fac = -vpdum*Dtemp*vdum     & 
                     & /(vdum**2 + vpdum**2)
                else
                  fac = 0.E0
                endif
                
                if(en_scatter)then
                  call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
                  fad = vpdum*vdum*Dtemp &
                      & /(vdum**2 + vpdum**2)
                else
                  fad = 0.E0
                end if
              end if
              
              elem%jloc=j+1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)       
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j+1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)   
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr) 
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

              elem%jloc=j-1
              elem%kloc=k-1
              call vgridboundary(2,elem,ingrida)
              if(ingrida)then
                jjh = indx(ifdis,imod,ix,i,elem%jloc,elem%kloc,is)
                mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
                elem%val=mat_elem
                call add_element(elem,ierr)
                if(neoclassics) then
                  call nc_copy_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)                
                end if
              end if

            end do
          end do
        end do

      end do
    end do
  end do

end subroutine collision_differential_numu

!--------------------------------------------------------------------
!> implements the simple model for momentum conservation
!> The integrals are performed in exp_integration to save time
!-------------------------------------------------------------------
subroutine cons_momentum

  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp
  use velocitygrid,   only : vpgr,mugr
  use dist,           only : fmaxwl,i_mom,i_ene,ifdis 
  use matdat,         only : add_element, set_indx
  use structures,     only : matrix_element
  use geom,           only : bn
  use specfun,        only : erf => sf_erf

  character(len=64) :: term='Collisions cons_mon'  
  character(len=64) :: term2='Collisions cons_ene'  

  integer :: i,j,k,ix,imod,is,ierr
  type (matrix_element) :: elem
  real :: vsqr
  real :: dumv

 
  if(mom_conservation)then
    ! Here is the momentum conservation term
    ! Identifier of the term 
    elem%term = term
    ! Type of iih and jjh 
    elem%itype = ifdis 
    elem%itloc = i_mom
    elem%ideriv = 0

    select case(cons_type)
      case('Lin')

       do is=1,nsp;do ix=1,nx;do imod=1,nmod

        do i=1,ns
         do j=1,nmu
          do k = 1,nvpar
           call set_indx(elem,imod,ix,i,j,k,is)
           vsqr = vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j)
           dumv = sqrt(vsqr)
           elem%val = fmaxwl(ix,i,j,k,is)*(erf(dumv)-dumv    &
           & *erfp(dumv))/dumv**3.*vpgr(i,j,k,is)
           call add_element(elem,ierr) 
          end do
         end do
       end do

      end do;end do;end do
     case('Xu')

      do is=1,nsp;do ix=1,nx;do imod=1,nmod

       do i=1,ns
        do j=1,nmu
         do k = 1,nvpar
          call set_indx(elem,imod,ix,i,j,k,is)
          elem%val = fmaxwl(ix,i,j,k,is)*vpgr(i,j,k,is)
          call add_element(elem,ierr) 
          end do
         end do
        end do

       end do;end do;end do
    end select
  end if

  if(ene_conservation)then
    ! Here is the energy conservation term
    ! Identifier of the term 
    elem%term = term2
    ! Type of iih and jjh 
    elem%itype = ifdis 
    elem%itloc = i_ene
    elem%ideriv = 0

    select case(cons_type)
      case('Lin')

       do is=1,nsp;do ix=1,nx;do imod=1,nmod

        do i=1,ns
         call maxwell_integrals(i,is)
         do j=1,nmu
          do k = 1,nvpar
           call set_indx(elem,imod,ix,i,j,k,is)
           vsqr = vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j)
           dumv = sqrt(vsqr)
           elem%val = fmaxwl(ix,i,j,k,is)*(erf(dumv)-2.*dumv  &
           & *erfp(dumv))/dumv
           call add_element(elem,ierr) 
          end do
         end do
        end do

       end do;end do;end do   
       case('Xu')

        do is=1,nsp;do ix=1,nx;do imod=1,nmod

         do i=1,ns
          call maxwell_integrals(i,is)
          do j=1,nmu
           do k = 1,nvpar
            call set_indx(elem,imod,ix,i,j,k,is)
            vsqr = vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j)
            elem%val = fmaxwl(ix,i,j,k,is)*(vsqr-(ene/partnum))
            call add_element(elem,ierr) 
           end do
          end do
         end do

        end do;end do;end do
    end select

  end if

end subroutine cons_momentum

!--------------------------------------------------------------------
!>Calculation of the prefactors of the integrals in the 
!>momentum conserving term
!-------------------------------------------------------------------
subroutine coll_mom_change_int

  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp,n_vpar_grid 
  use grid,           only : gvpar, gmu, n_mu_grid
  use velocitygrid,   only : mugr,vpgr,dvp,dmu,vgridboundaryMom
  use dist,           only : ifdis
  use index_function, only : indx
  use geom,           only : bn 
  use components,     only : vthrat
  use general,        only : gkw_abort 
  use control,        only : flux_tube

  character(len=64) :: term='Collisions coll_mom_change_int'

  integer :: i,j,k,ix,imod,is
  integer :: kdum,jdum
  real    :: Dtemp, vdum, mudum,vsqr
  logical :: ingrida
  real    :: mat_elem
  integer :: jjh
  real    :: fac, fad, faf

  if (.not. flux_tube) then
    ! fmaxwl is calculated at the first nx grid point  
    call gkw_abort('Collisions conservation has not been made global - yet.')
  endif 

  do i=1,ns; do is=1,nsp

    call maxwell_integrals(i,is)

    do ix=1,nx ; do imod = 1,nmod 

      do j=0,nmu+1 ; do k=0,nvpar+1
        
        vsqr = vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j)
        
        !First is the dvpardvpar term            
        vdum = vpgr(i,j,k,is)+0.5E0*dvp
        mudum = mugr(j)
        if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
          fac = 0.E0
          fad = 0.E0
          faf = 0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            fac = 2.E0*mudum*Dtemp*bn(ix,i)       & 
              & /(vdum**2 + 2.E0*mudum*bn(ix,i))
            fac = fac/(dvp*dvp)
          else
            fac=0.E0
          endif
                
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fad = vdum**2*Dtemp/sqrt(vdum**2+2*mudum*bn(ix,i)) 
            fad = fad/(dvp*dvp)                  
          else
            fad=0.E0
          endif
                
          if(friction_coll)then
            if(selfcollcon)then
              call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            faf = vdum*Dtemp/(dvp)
          else
            faf=0.E0
          endif
        endif
    
        jdum=j
        kdum=k+1
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(fac+fad+0.5E0*faf)  
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(-fac-fad+0.5E0*faf) 
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        vdum = vpgr(i,j,k,is)-0.5*dvp
        mudum = mugr(j)     
        if((gvpar(k).eq.1).and.mass_conserve)then
          fac = 0.E0
          fad = 0.E0
          faf = 0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = 2.E0*mudum*Dtemp*bn(ix,i)       & 
                  & /(vdum**2 + 2*mudum*bn(ix,i))
            fac = fac/(dvp*dvp)
          else
            fac = 0.E0
          endif
                
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fad = vdum**2*Dtemp/sqrt(vdum**2+2*mudum*bn(ix,i)) 
            fad = fad/(dvp*dvp)     
          else
            fad = 0.E0
          endif
                
          if(friction_coll)then
            if(selfcollcon)then
              call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            faf = vdum*Dtemp/(dvp)     
          else
            faf = 0.E0
          endif
        end if
  
        jdum=j
        kdum=k
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(-fac-fad-0.5E0*faf) 
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k-1
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(fac+fad-0.5E0*faf) 
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        !Second is the dmudmu term    
   
        vdum = vpgr(i,j,k,is)
        mudum= mugr(j)+0.5E0*dmu
        if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
          fac = 0.E0
          fad = 0.E0
          faf = 0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)    
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = 2.E0*mudum*Dtemp*vdum**2     & 
                  & /((vdum**2 + 2.E0*mudum*bn(ix,i))*bn(ix,i))
            fac = fac/(dmu*dmu)
          else
            fac = 0.E0
          endif
                
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            fad = 4.E0*mudum**2*Dtemp/sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
            fad = fad/(dmu*dmu)
          else
            fad = 0.E0
          end if
                  
          if(friction_coll)then
            if(selfcollcon)then
              call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            faf = 2.E0*mudum*Dtemp/(dmu)
          else
            faf = 0.E0
          endif
        endif
          
        jdum=j+1
        kdum=k               
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida) 
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(fac+fad+0.5E0*faf)  
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(-fac-fad+0.5E0*faf)  
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        vdum = vpgr(i,j,k,is)
        mudum = mugr(j)-0.5E0*dmu
        if((gmu(j).eq.1).and.mass_conserve)then
          fac = 0.E0
          fad = 0.E0
          faf = 0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = 2.E0*mudum*Dtemp*vdum**2     &
                  & /((vdum**2 + 2.E0*mudum*bn(ix,i))*bn(ix,i))
            fac = fac/(dmu*dmu)
          else
            fac = 0.E0
          endif  
                
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            fad = 4.E0*mudum**2*Dtemp/sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
            fad = fad/(dmu*dmu)
          else
            fad = 0.E0
          endif
                
          if(friction_coll)then
            if(selfcollcon)then
              call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            faf = 2.E0*mudum*Dtemp/(dmu)
          else
            faf = 0.E0
          endif
        endif

        jdum=j
        kdum=k
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then               
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(-fac-fad-0.5E0*faf)  
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k   
        call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = vthrat(is)*(fac+fad-0.5E0*faf) 
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if
          
        !Cross terms (there are no cross terms in friction)
        !Must interpolate (using Bilinear) the half/half grid
        !points

        !Firstly dmudvpar
        vdum = vpgr(i,j,k,is)
        mudum = mugr(j)+0.5*dmu
        if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
          fac=0.E0
          fad=0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = -2*mudum*Dtemp*vdum     & 
                    & /(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fac = 0.E0
          endif
                 
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            fad = 2.E0*mudum*vdum*Dtemp &
                &  /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fad = 0.E0
          end if
        endif
          
        jdum=j+1
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)                  
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j+1
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)                     
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        vdum = vpgr(i,j,k,is)
        mudum = mugr(j)-0.5E0*dmu        
        if((gmu(j).eq.1).and.mass_conserve)then
          fac=0.E0
          fad=0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = -2*mudum*vdum*Dtemp     & 
                & /(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fac = 0.E0
          endif
                 
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fad = 2.E0*mudum*vdum*Dtemp &
                &  /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fad = 0.E0
          end if
        end if
        
        jdum=j
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)  
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh =indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        !Then dvpardmu
        vdum = vpgr(i,j,k,is)+0.5E0*dvp
        mudum = mugr(j)
        if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
          fac=0.E0
          fad=0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = -2.E0*mudum*Dtemp*vdum     & 
                    & /(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fac = 0.E0
          endif
                
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            endif
            fad = 2.E0*mudum*vdum*Dtemp &
                  & /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fad = 0.E0
          end if
        end if

        jdum=j+1
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j+1
        kdum=k
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k+1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)     
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)    
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        vdum = vpgr(i,j,k,is)-0.5*dvp
        mudum = mugr(j)
        if((gvpar(k).eq.1).and.mass_conserve)then
          fac=0.E0
          fad=0.E0
        else
          if(pitch_angle)then
            if(selfcollcon)then
              call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            else
              call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fac = -2*mudum*Dtemp*vdum     & 
                  & /(vdum**2 + 2*mudum*bn(ix,i))
          else
            fac = 0.E0
          endif
               
          if(en_scatter)then
            if(selfcollcon)then
              call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
            else
              call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
            endif
            fad = 2.E0*mudum*vdum*Dtemp &
               & /sqrt(vdum**2 + 2.E0*mudum*bn(ix,i))
          else
            fad = 0.E0
          end if
        end if
         
        jdum=j+1
        kdum=k
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)       
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j+1
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)   
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

        jdum=j-1
        kdum=k-1
        call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
        if(ingrida)then
          jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
          mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dmu)
          call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
        end if

      end do; end do !vpar,mu

    end do; end do   !ix,imod 

  end do; end do      !is,i

  !if (allocated(gammab)) deallocate(gammab)

end subroutine coll_mom_change_int

!--------------------------------------------------------------------
!>Calculation of the prefactors of the integrals in the 
!>momentum conserving term
!-------------------------------------------------------------------
subroutine coll_mom_change_int_numu
  
  use grid,           only : nmod,nx,ns,nmu,nvpar,nsp,gvpar
  use grid,           only : n_vpar_grid, gmu, n_mu_grid
  use velocitygrid,   only : mugr,vpgr, dvp, dvperp, vgridboundaryMom
  use dist,           only : ifdis
  use index_function, only : indx
  use geom,           only : bn 
  use components,     only : vthrat
  
  character(len=64) :: term='Collisions coll_mom_change_int_numu'

  integer :: i,j,k,ix,imod,is
  integer :: kdum,jdum
  real    :: Dtemp, vdum, mudum, vpdum
  real    :: dvrp
  logical :: ingrida
  real    :: mat_elem
  integer :: jjh
  real    :: fac, fad, faf, vzero, del, vsqr

  !if (non_local_profiles) then
  !  ! fmaxwl is calculated at the first nx grid point  
  !  write(*,*)'Momentum conservation has not been made global - yet!'
  !  call gkw_abort('Sorry')
  !endif 

  do i=1,ns; do is=1,nsp

    call maxwell_integrals(i,is)
    
    do ix=1,nx ; do imod = 1,nmod 

      do j=0,nmu+1 ; do k=0,nvpar+1

         vsqr = vpgr(i,j,k,is)**2 + 2.E0*bn(ix,i)*mugr(j)
         vzero = sqrt(2.E0*mugr(j)*bn(ix,i))
         dvrp = sqrt(bn(ix,i))*dvperp
       
         !First is the dvpardvpar term            
         vdum = vpgr(i,j,k,is)+0.5E0*dvp
         mudum = mugr(j)
         vpdum = sqrt(2.E0*bn(ix,i)*mudum)
         if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
           fac = 0.E0
           fad = 0.E0
           faf = 0.E0
         else
           if(pitch_angle)then
             if(selfcollcon)then
               call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             fac = vpdum**2*Dtemp       & 
                  & /(vdum**2 + vpdum**2)
             fac = fac/(dvp*dvp)
           else
             fac=0.E0
           endif
               
           if(en_scatter)then
             if(selfcollcon)then
               call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             fad = vdum**2*Dtemp/(vdum**2+vpdum**2)
             fad = fad/(dvp*dvp)                  
           else
             fad=0.E0
           endif
                
           if(friction_coll)then
             if(selfcollcon)then
               call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             faf = vdum*Dtemp/(sqrt(vdum**2+vpdum**2)*dvp)
           else
             faf=0.E0
           endif

         endif

         call ccdelta(fac+fad,faf,dvp,del)
              
         jdum=j
         kdum=k+1
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(fac+fad+(1.E0-del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if

         jdum=j
         kdum=k
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(-fac-fad+del*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if
     
         vdum = vpgr(i,j,k,is)-0.5*dvp
         mudum = mugr(j)    
         vpdum = sqrt(2.E0*bn(ix,i)*mudum)
         if((gvpar(k).eq.1).and.mass_conserve)then
           fac = 0.E0
           fad = 0.E0
           faf = 0.E0
         else
           if(pitch_angle)then
             if(selfcollcon)then
               call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             fac = vpdum**2*Dtemp       & 
                & /(vdum**2 + vpdum**2)
             fac = fac/(dvp*dvp)
           else
             fac = 0.E0
           endif
                
           if(en_scatter)then
             if(selfcollcon)then
               call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             fad = vdum**2*Dtemp/(vdum**2+vpdum**2) 
             fad = fad/(dvp*dvp)     
           else
              fad = 0.E0
           endif
                
           if(friction_coll)then
             if(selfcollcon)then
               call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             faf = vdum*Dtemp/(sqrt(vdum**2+vpdum**2)*dvp)     
           else
             faf = 0.E0
           endif

         end if
                 
         jdum=j
         kdum=k
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if
         
         jdum=j
         kdum=k-1
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if

         !Second is the dvperpdvperp term    
           
         vdum = vpgr(i,j,k,is)
         vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) + 0.5E0*dvrp
         mudum = vpdum**2/(2*bn(ix,i))
              
         if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
           fac = 0.E0
           fad = 0.E0
           faf = 0.E0
         else
           if(pitch_angle)then
             if(selfcollcon)then
               call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)    
             else
               call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)    
             endif
             fac = vpdum*Dtemp*vdum**2     & 
                & /(vdum**2 + vpdum**2)
             fac = fac/(vzero*dvrp*dvrp)
           else
             fac = 0.E0
           endif
                
           if(en_scatter)then
             if(selfcollcon)then
               call selfcaldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i) 
             else
               call caldvv(vpgr(i,j,k,is),mudum*bn(ix,i),ix,is,Dtemp,i)             
             endif
             fad = vpdum**3*Dtemp/(vdum**2 + vpdum**2)
             fad = fad/(vzero*dvrp*dvrp)
           else
             fad = 0.E0
           end if
                  
           if(friction_coll)then
             if(selfcollcon)then
               call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             else
               call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
             endif
             faf = vpdum**2*Dtemp/(vzero*sqrt(vdum**2 + vpdum**2)*dvrp)
           else
             faf = 0.E0
           endif

         endif
         
         call ccdelta(fac+fad,faf,dvrp,del)
      
         jdum=j+1
         kdum=k               
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida) 
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(fac+fad+ (1.E0-del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if

         jdum=j
         kdum=k
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(-fac-fad+del*faf) 
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if
      
         vdum = vpgr(i,j,k,is)
         vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) - 0.5E0*dvrp
         mudum = vpdum**2/(2*bn(ix,i))
              
         !No need for imposition of mass-conservation on the lower
         !mu boundary as the flux is naturally zero at this point.
         if(pitch_angle)then
           if(selfcollcon)then
             call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           else
             call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           endif
           fac = vpdum*Dtemp*vdum**2     &
              & /(vdum**2 + vpdum**2)
           fac = fac/(vzero*dvrp*dvrp)
         else
           fac = 0.E0
         endif  
                
         if(en_scatter)then
           if(selfcollcon)then
             call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           else
             call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           endif
           fad = vpdum**3*Dtemp/(vdum**2 + vpdum**2)
           fad = fad/(vzero*dvrp*dvrp)
         else
           fad = 0.E0
         endif
                
         if(friction_coll)then
           if(selfcollcon)then
             call selfcalfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           else
             call calfv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
           endif
           faf = vpdum**2*Dtemp/(dvrp*sqrt(vdum**2 + vpdum**2)*vzero)
         else
           faf = 0.E0
         endif

         call ccdelta(fac+fad,faf,dvrp,del)

         jdum=j
         kdum=k
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then               
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(-fac-fad-(del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if

         jdum=j-1
         kdum=k   
         call vgridboundaryMom(1,jdum,kdum,j,k,ingrida)
         if(ingrida)then
           jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
           mat_elem = vthrat(is)*(fac+fad-(1.E0-del)*faf)
           call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
         end if
     
         !Cross terms (there are no cross terms in friction)
         !Must interpolate (using Bilinear) the half/half grid
         !points        

          !Firstly dvperpdvpar
          vdum = vpgr(i,j,k,is)
          vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) + 0.5E0*dvrp
          mudum = vpdum**2/(2.E0*bn(ix,i))
          if((gmu(j).eq.n_mu_grid).and.mass_conserve)then
            fac=0.E0
            fad=0.E0
          else
            if(pitch_angle)then
              if(selfcollcon)then
                call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              else
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              endif
              fac = -vpdum**2*Dtemp*vdum     & 
               & /(vzero*(vdum**2 + vpdum**2))
            else
              fac = 0.E0
            endif
                 
            if(en_scatter)then
              if(selfcollcon)then
                call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              else
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              endif
              fad = vpdum**2*vdum*Dtemp &
                &  /(vzero*(vdum**2 + vpdum**2))
            else
              fad = 0.E0
            end if

          end if
              
          jdum=j+1
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)  
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j+1
          kdum=k-1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j
          kdum=k-1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if


          vdum = vpgr(i,j,k,is)
          vpdum = sqrt(2.E0*bn(ix,i)*mugr(j)) - 0.5E0*dvrp
          mudum = vpdum**2/(2*bn(ix,i))    
          if((gmu(j).eq.1).and.mass_conserve)then
            fac=0.E0
            fad=0.E0
          else
            if(pitch_angle)then
              if(selfcollcon)then
                call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              else
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              endif
              fac = -vpdum**2*vdum*Dtemp     & 
                  & /(vzero*(vdum**2 + vpdum**2))
            else
              fac = 0.E0
            endif
                 
            if(en_scatter)then
              if(selfcollcon)then
                call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              else
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              endif
              fad = vpdum**2*vdum*Dtemp &
                 &  /(vzero*(vdum**2 + vpdum**2))
            else
              fad = 0.E0
            end if
          end if
            
          jdum=j
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j-1
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j-1
          kdum=k-1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida) 
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j
          kdum=k-1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh =indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          !Then dvpardvperp
          vdum = vpgr(i,j,k,is)+0.5E0*dvp
          mudum = mugr(j)
          vpdum = sqrt(2.E0*bn(ix,i)*mugr(j))
          if((gvpar(k).eq.n_vpar_grid).and.mass_conserve)then
            fac=0.E0
            fad=0.E0
          else
            if(pitch_angle)then
              if(selfcollcon)then
                call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              else
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              endif
              fac = -vpdum*vdum*Dtemp     & 
                 & /(vdum**2 + vpdum**2)
            else
               fac = 0.E0
            endif
                 
            if(en_scatter)then
              if(selfcollcon)then
                call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              else
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              endif
              fad = vpdum*vdum*Dtemp &
                 & /(vdum**2 + vpdum**2)
            else
              fad = 0.E0
            end if
          end if
    
          jdum=j+1
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j+1
          kdum=k
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j-1
          kdum=k+1
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)     
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          jdum=j-1
          kdum=k
          call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)   
          if(ingrida)then
            jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
            mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
            call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
          end if

          vdum = vpgr(i,j,k,is)-0.5*dvp
          mudum = mugr(j)
          vpdum = sqrt(2.E0*bn(ix,i)*mugr(j))
          if((gvpar(k).eq.1).and.mass_conserve)then
            fac=0.E0
            fad=0.E0
          else
            if(pitch_angle)then
              if(selfcollcon)then
                call selfcaldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              else
                call caldthth(vdum,mudum*bn(ix,i),ix,is,Dtemp,i)
              endif
              fac = -vpdum*Dtemp*vdum     & 
                  & /(vdum**2 + vpdum**2)
            else
              fac = 0.E0
            endif
                
            if(en_scatter)then
              if(selfcollcon)then
                call selfcaldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              else
                call caldvv(vdum,mudum*bn(ix,i),ix,is,Dtemp,i) 
              endif
              fad = vpdum*vdum*Dtemp &
                    & /(vdum**2 + vpdum**2)
              else
                fad = 0.E0
              end if
            end if
              
            jdum=j+1
            kdum=k
            call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)       
            if(ingrida)then
              jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
              mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
              call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
            end if

            jdum=j+1
            kdum=k-1
            call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)   
            if(ingrida)then
              jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
              mat_elem = -0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
              call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
            end if

            jdum=j-1
            kdum=k
            call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
            if(ingrida)then
              jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
              mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
              call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
            end if

            jdum=j-1
            kdum=k-1
            call vgridboundaryMom(2,jdum,kdum,j,k,ingrida)
            if(ingrida)then
              jjh = indx(ifdis,imod,ix,i,jdum,kdum,is)
              mat_elem = 0.25*vthrat(is)*(fac+fad)/(dvp*dvrp)
              call conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)
            end if

      end do; end do !vpar,mu

    end do; end do   !ix,imod 

  end do; end do      !is,i
  
  !if (allocated(gammab)) deallocate(gammab)
  
end subroutine coll_mom_change_int_numu


!****************************************************************************
!> Derivative of erf  \f[ d (erf)/dx = 2*exp(-x^2)/ sqrt(\pi) \f]
!----------------------------------------------------------------------------
function erfp(x)

  real, intent(in) :: x
  real :: erfp
  real, parameter :: two_over_rootpi = 1.12837916709551257389615890E0
  real, parameter :: rootln10 = 1.51742712938514635086297239E0

  if (x > rootln10*sqrt(0.99*range(x))) then
    erfp = 0. ! don't bother calculating for x large.
  else
    erfp = two_over_rootpi*exp(-x**2)
  end if

end function erfp

!****************************************************************************
!Calculates delta, the Chang-Cooper differencing parameter for the differential
!terms in the collision operator.
!----------------------------------------------------------------------------
subroutine ccdelta(fac,faf,dv,delta)

  real, intent(in) :: fac,faf,dv
  real, intent(out) :: delta
  real :: dum
 
  !At the moment set to 0.5 always (i.e centrally differenced.
  delta = 0.5E0
  return

  if(abs(fac).lt.(1.0E-10))then
    delta = 0.E0
  else if(abs(faf).lt.(1.0E-10))then
    delta = 0.5E0
  else
    dum = -faf*dv/fac
    delta = (1.E0/dum) - (1.E0/(exp(dum)-1.E0))
  end if
  
end subroutine ccdelta

!****************************************************************************
!> Copy the collision operator matrix entries to the appropriate matrices for
!> the neoclassical diagnostis
!----------------------------------------------------------------------------
subroutine nc_copy_matelem (term,mat_elem_nc,jjh, imod,ix, i, j, k,is  )
  use matrix_format, only : put_element
  use geom,       only : bn, ints, efun, signB, bt_frac, Rfun, jfun
  use components, only : tmp, mas, signz
  use matdat,     only : matn,matn_e,matn_v
  use velocitygrid, only : mugr, vpgr, intvp, intmu
  use mode,     only : ixzero, iyzero 
  use constants, only : pi
  use rotation,       only : cf_trap, cf_drift, cfen, vcor
  !use linear_terms,   only : lneorotsource
  integer, intent(in) :: jjh, imod, ix, i, j, k, is
  complex, intent(in) :: mat_elem_nc
  character(len=64), intent(in) :: term
  complex :: dum_mat_elem, dum_mat_elem2
  real :: dum, velshift

  if ((ix .eq. ixzero) .and. (imod .eq. iyzero)) then
        
        velshift = signB*bt_frac(ix,i) * rfun(ix,i)*vpgr(i,j,k,is)
        !if(.not.lneorotsource)then
          velshift = velshift + vcor*jfun(ix,i)
        !endif
        dum_mat_elem = mat_elem_nc * velshift*4*pi*efun(ix,i,1,2)*bn(ix,i)/signz(is)* &
                   & sqrt(mas(is)*tmp(ix,is))*intmu(j)* intvp(i,j,k,is)*ints(i)
        dum = (vpgr(i,j,k,is)**2 +  2.E0*mugr(j)*bn(ix,i))
        !Correction for when the centrifugal drift is included.
        if(cf_drift.or.cf_trap)then
          dum = dum + cfen(i,is)
        endif

        call put_element(matn,1,jjh,dum_mat_elem)
        dum_mat_elem2 = dum_mat_elem * (dum-5.0E0/2.0E0)
        call put_element(matn_e,1,jjh,dum_mat_elem2)
        dum_mat_elem2 = dum_mat_elem *0.5E0*velshift
        call put_element(matn_v,1,jjh,dum_mat_elem2)
  end if

end subroutine nc_copy_matelem

!****************************************************************************
!> Places the integral conservation terms in the matrix
!----------------------------------------------------------------------------
subroutine conserve_put_matelem(term,mat_elem,jjh,imod,ix,i,j,k,is)

  use geom,       only : bn
  use matdat,     only : matm
  use matrix_format, only : put_element
  use velocitygrid, only : mugr, vpgr, intvp, intmu
  use dist,       only : i_mom,i_ene
  use index_function, only : indx

  integer, intent(in) :: jjh, imod, ix, i, j, k, is
  real, intent(in) :: mat_elem
  !> FIXME the matrix is real values, but at the moment, there is only
  !> a complex sparse matrix datatype
  !real :: mat_elem_m
  complex :: mat_elem_m
  real :: vsqr
  integer :: iih
  character(len=64), intent(in) ::term

  if(mom_conservation)then
    select case(cons_type)
     case('Lin')
      iih = indx(i_mom,imod,ix,i,is)
      mat_elem_m = -mat_elem*vpgr(i,j,k,is)*intvp(i,j,k,is)*intmu(j)/erfintm
      call put_element(matm,iih,jjh,mat_elem_m)
     case('Xu')
      iih = indx(i_mom,imod,ix,i,is)
      mat_elem_m = -mat_elem*vpgr(i,j,k,is)*intvp(i,j,k,is)*intmu(j)/parmom
      call put_element(matm,iih,jjh,mat_elem_m)
    end select
  end if

  if(ene_conservation)then
    select case(cons_type)
     case('Lin')
      iih = indx(i_ene,imod,ix,i,is)
      vsqr = vpgr(i,j,k,is)**2 + 2*bn(ix,i)*mugr(j)
      mat_elem_m = -mat_elem*vsqr*intvp(i,j,k,is)*intmu(j)/erfinte     
      call put_element(matm,iih,jjh,mat_elem_m)
     case('Xu')
      iih = indx(i_ene,imod,ix,i,is)
      vsqr = vpgr(i,j,k,is)**2 + 2*bn(ix,i)*mugr(j)
      mat_elem_m = -mat_elem*vsqr*intvp(i,j,k,is)*intmu(j)/(enesqr-(ene*ene/partnum)) 
      call put_element(matm,iih,jjh,mat_elem_m)
    end select
  end if 

end subroutine conserve_put_matelem

!****************************************************************************
!> Routine that performs all the integrals over the Maxwellian needed for the
!> conservation terms
!----------------------------------------------------------------------------
subroutine maxwell_integrals(i,is) 

  use geom,         only : bn
  use velocitygrid, only : mugr, vpgr, intvp, intmu
  use dist,         only : fmaxwl
  use grid,         only : nvpar, nmu
  use mpiinterface, only : mpiallreduce_sum
  use mpicomms,     only : COMM_VPAR_NE_MU_NE
  use specfun,      only : erf => sf_erf

  real :: fac 
  integer, intent(in) :: i,is
  integer :: j,k
  real :: dum,dume,dum4,dumn
  real :: dumerfm, dumerfe, dumv

  ! For each i and is, calculate the integral over the maxwellian first.
  dum  = 0.E0
  dume = 0.E0
  dum4 = 0.E0
  dumn = 0.E0
  dumerfm = 0.E0
  dumerfe = 0.E0

  do j=1,nmu
    do k=1,nvpar
      dum = dum + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)*vpgr(i,j,k,is)**2   
      dumn = dumn + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)
      fac = vpgr(i,j,k,is)**2 + 2.E0*bn(1,i)*mugr(j)
      dumv = sqrt(fac)
      dum4 = dum4 + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)*fac*fac
      dume = dume + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)*(fac)   
      dumerfm = dumerfm + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)*          &
      & (erf(dumv)-dumv*erfp(dumv))/dumv**3.*vpgr(i,j,k,is)**2.
      dumerfe = dumerfe + intvp(i,j,k,is)*intmu(j)*fmaxwl(1,i,j,k,is)*          &
      & (erf(dumv)-2.*dumv*erfp(dumv))/dumv*fac
    end do
  end do

  ! reduce over points with the same every apart from vpar and mu
  ! WARNING dum*2 are global variables
  call mpiallreduce_sum(dum,parmom,1,COMM_VPAR_NE_MU_NE)
  call mpiallreduce_sum(dume,ene,1,COMM_VPAR_NE_MU_NE)
  call mpiallreduce_sum(dumn,partnum,1,COMM_VPAR_NE_MU_NE)
  call mpiallreduce_sum(dum4,enesqr,1,COMM_VPAR_NE_MU_NE)
  call mpiallreduce_sum(dumerfe,erfinte,1,COMM_VPAR_NE_MU_NE)
  call mpiallreduce_sum(dumerfm,erfintm,1,COMM_VPAR_NE_MU_NE)

end subroutine maxwell_integrals

end module collisionop
