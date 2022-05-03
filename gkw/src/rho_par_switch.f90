!-------------------------------------------------------------------------------
!> Contains switches and namelist for parallel derivative rho star effects
!-------------------------------------------------------------------------------
module rho_par_switch 

  implicit none 

  public :: rho_par_switch_write_nml, rho_par_switch_read_nml 
  public :: rho_par_switch_bcast_nml, rho_par_switch_check_params
  
  public :: ltrapdf_rhostar, lvdgrad_phi_fm_rhostar, s_average 
  public :: lve_grad_fm_rhostar, lvdgradf_rhostar, lnonlinear_rhostar
  public :: lflux_rhostar 

  private 

  !
  !The switches of the finite rho* parallel derivative terms.
  !

  !> Logical the switches on/off the rho* correction in the trapping term  
  logical, save :: ltrapdf_rhostar        = .false. 

  !> Logical that switches on/off the rho* correction on the drift in the 
  !> gradient of the potential 
  logical, save :: lvdgrad_phi_fm_rhostar = .false.

  !> Logical that switches on/off the rho* correction on the ExB  drift in 
  !> the gradient of the Maxwell 
  logical, save :: lve_grad_fm_rhostar    = .false.

  !> Logical that switches on/off the rho* correction on the drift in the 
  !> gradient of the perturbed distribution. 
  logical, save :: lvdgradf_rhostar       = .false.

  !> Logical to enable the finite rhostar direction in the gyro average and 
  !> polarisation. This option exclueds the useage of parallel s and Hermitian 
  !> operator.
  logical, save :: s_average              = .false. 

  !> logical that determines if the finite rho* parallel derivatives are taken 
  !> into account in the non-linear terms. 
  logical, save :: lnonlinear_rhostar     = .false. 
 
  !> logical that deterimes if the finite rho* parallel derivative is taken 
  !> into account in the calculation of the fluxes 
  logical, save :: lflux_rhostar          = .false. 

  interface rho_par_switch_write_nml
     module procedure rho_par_switch_read_nml
  end interface

contains 

subroutine rho_par_switch_read_nml(lun,io_stat,lwrite)
  use io, only : write_run_parameter
  use mpiinterface, only : root_processor
  integer, intent(in)  :: lun
  integer, intent(out) :: io_stat

  logical, optional, intent(in) :: lwrite

  namelist /finite_rho_parallel / ltrapdf_rhostar, lvdgrad_phi_fm_rhostar, &
           & lve_grad_fm_rhostar, lvdgradf_rhostar, lnonlinear_rhostar,    & 
           & s_average, lflux_rhostar
          
  io_stat = 0
  if (present(lwrite)) then
    if (.not. lwrite) then       
      read(lun,NML=finite_rho_parallel,IOSTAT=io_stat) 
    else
      ! do nothing
    end if
  else
    if(root_processor) write(lun,NML=finite_rho_parallel)

    call write_run_parameter('finite_rho_parallel', 'ltrapdf_rhostar', ltrapdf_rhostar)
    call write_run_parameter('finite_rho_parallel', 'lvdgrad_phi_fm_rhostar', lvdgrad_phi_fm_rhostar)
    call write_run_parameter('finite_rho_parallel', 'lve_grad_fm_rhostar', lve_grad_fm_rhostar)
    call write_run_parameter('finite_rho_parallel', 'lvdgradf_rhostar', lvdgradf_rhostar)
    call write_run_parameter('finite_rho_parallel', 'lnonlinear_rhostar', lnonlinear_rhostar)
    call write_run_parameter('finite_rho_parallel', 's_average', s_average)
    call write_run_parameter('finite_rho_parallel', 'lflux_rhostar', lflux_rhostar)

  end if 

end subroutine rho_par_switch_read_nml

!------------------------------------------------------------------------------
!> check the switches namelist params
!------------------------------------------------------------------------------
subroutine rho_par_switch_check_params

  use grid,       only : n_procs_s
  use general,    only : gkw_exit
  
  if(lflux_rhostar.and.(n_procs_s.gt.1))then
    call gkw_exit('There is currently a bug using the finite rhostar with parallel s')
  endif

end subroutine rho_par_switch_check_params

!------------------------------------------------------------------------------
!> bcast the switches namelist params
!------------------------------------------------------------------------------
subroutine rho_par_switch_bcast_nml

  use mpiinterface, only : mpibcast 
  
  call mpibcast(ltrapdf_rhostar,        1)
  call mpibcast(lvdgrad_phi_fm_rhostar, 1)
  call mpibcast(lve_grad_fm_rhostar,    1)
  call mpibcast(lvdgradf_rhostar,       1)
  call mpibcast(s_average           ,   1) 
  call mpibcast(lnonlinear_rhostar,     1) 
  call mpibcast(lflux_rhostar,          1)

end subroutine rho_par_switch_bcast_nml

end module rho_par_switch 
