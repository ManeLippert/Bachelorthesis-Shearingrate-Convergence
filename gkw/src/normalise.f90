!-----------------------------------------------------------------------------
!> This module contains a routine to re-normalise the solution in linear runs
!> so that arithmetic overflows do not occur.
!>
!> The amplitudes used to determine the normalisation factors are public so that
!> diagnostics can access them, e.g. to compute growth rates.
!>
!> This module also has the effect of normalising linear fluxes to |phi|^2
!> which is useful for quasilinear studies.
!>
!-----------------------------------------------------------------------------
module normalise

  implicit none

  private

  public :: calc_norm_factor
  public :: normalise_fdisi, normalise_fdisk_with_fdisi, normalise_init
  public :: rotate_fdis, get_cmplx_rotation_factor, normalise_after_timestep

  !> amplitude
  real, save, public :: amp
  !> total accumulated single normalisation factor
  real, allocatable, save, public :: accumulated_normfactor(:)

  !> normalisation factor for each toroidal mode
  real, save, public, allocatable :: amp_per_mode(:), fnorm1d(:)

  !> ready to use switch
  logical, save :: initialised = .false.

  ! arrays for index
  integer, save, allocatable :: findx(:,:), lindx(:,:)
  
  ! scaling of the normalised distribution function for nonlinear runs
  real, save, public :: nonlin_norm_fac

contains


!-------------------------------------------------------------------------------
!> Setup everything necessary to use the other routines in the module.
!> - Allocate the necessary arrays
!> - Create an index array to reference the fields for each mode. For a
!>   single normalisation factor this is not necessary as the fields can be
!>   accessed as one block.
!> - Create an index array to reference each mode of fdisi to be normalised.
!>   With a single factor, all fdisi is normalised with that factor.
!------------------------------------------------------------------------------

subroutine normalise_init()
  use mode, only : ikxspace
  use control,        only : nlphi, nlapar, nlbpar, normalized, non_linear
  use control,        only : normalize_per_toroidal_mode, spectral_radius
  use grid,           only : nx, nmod, ns, nsp, nmu, nvpar
  use dist,           only : ifdis, iphi, iapar, ibpar, iphi_ga, iapar_ga, ibpar_ga, i_mom, i_ene, nf, nsolc
  use index_function, only : indx
  use general,        only : gkw_warn, gkw_abort
  use index_function, only : IS_3D_FIELD, IS_COLL_CONS_FIELD, IS_GYROAVG_FIELD
  use collisionop, only : mom_conservation
  
  integer, allocatable :: field_id(:)
  integer :: ierr, idx, ipar, is, jv, kt, ix, imod, k, f_id
  integer :: number_of_all_fields

  allocate(fnorm1d(nmod),stat=ierr)
  if (ierr /= 0) call gkw_abort('normalise_init :: fnorm1d')
  allocate(accumulated_normfactor(nmod),stat=ierr)
  if (ierr /= 0) call gkw_abort('normalise_init :: accumulated_normfactor')
    
  if (normalize_per_toroidal_mode) then

    if (non_linear) then
      call gkw_warn('Control: Run is nonlinear but normalised&
         & (per toroidal mode)! This is not tested!')
      !> For nonlinear runs, saturation occurs and normalisation is not needed.
      !> But for special NL studies it is possible to normalise per toroidal mode.
      if(ikxspace > 1) then
        call gkw_warn('control: ikxspace should be 1 to use &
           & normalize_per_toroidal_mode in a non_linear run.')
      end if
    end if
  
    ! allocate an array for the normalisation factors for all toroidal
    allocate(amp_per_mode(nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('normalise_init :: amp_per_mode')
    ! allocate the index array for normalising
    allocate(findx(nmod,nf/nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('normalise_init :: findx')

    ! fill the index array
    idx = 0
    do is=1,nsp ; do kt=1,nvpar ; do jv=1,nmu ; do ipar=1,ns ; do ix=1,nx
      idx = idx + 1
      do imod = 1, nmod
        findx(imod,idx) = indx(ifdis,imod,ix,ipar,jv,kt,is)
      end do
    end do ; end do ; end do ; end do ; end do

    ! Construct the array for quick look up of the fields when calculating the
    ! normalisation factor.

    ! at max, we have 3 regular fields, 3 gyroavg fields and 2
    ! collisionop-related 'fields'
    allocate(field_id(3*2+2),stat=ierr)
    if (ierr /= 0) call gkw_abort('normalise_init :: field_id')
    number_of_all_fields = 0
    if (nlphi) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = iphi
    end if
    if (nlapar) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = iapar
    end if
    if (nlbpar) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = ibpar
    end if
    
    if(nlphi.and.(.not. spectral_radius)) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = iphi_ga
    end if
    if(nlapar.and.(.not. spectral_radius)) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = iapar_ga
    end if
    if(nlbpar.and.(.not. spectral_radius)) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = ibpar_ga
    end if
    if(mom_conservation) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = i_mom
    end if
    if(mom_conservation) then
      number_of_all_fields = number_of_all_fields + 1
      field_id(number_of_all_fields) = i_ene
    end if

    ! allocate the array
    allocate(lindx(nmod,(nsolc-nf)/nmod))
    if (ierr /= 0) call gkw_abort('normalise_init :: lindx')
    ! fill the array
    idx = 0
    lindx(:,:) = 0
    do k = 1, number_of_all_fields
      f_id = field_id(k)
      if(iand(f_id,IS_3D_FIELD)/=0) then
        do ipar = 1, ns
          do ix = 1, nx
            ! N.B. use nmod for inner loop, both here and elsewhere in the module
            idx = idx + 1
            do imod = 1, nmod
              lindx(imod,idx) = indx(f_id,imod,ix,ipar)
            end do
          end do
        end do
      else if(iand(f_id,IS_GYROAVG_FIELD)/=0) then
        do is = 1, nsp
          do jv = 1, nmu
            do ipar = 1, ns
              do ix = 1, nx
                ! N.B. use nmod for inner loop, both here and
                ! elsewhere in the module
                idx = idx + 1
                do imod = 1, nmod
                  lindx(imod,idx) = indx(f_id,imod,ix,ipar,jv,is)
                end do
              end do
            end do
          end do
        end do
      else if(iand(f_id,IS_COLL_CONS_FIELD)/=0) then
        do is = 1, nsp
          do ipar = 1, ns
            do ix = 1, nx
              ! N.B. use nmod for inner loop, both here and elsewhere
              ! in the module
              idx = idx + 1
              do imod = 1, nmod
                lindx(imod,idx) = indx(f_id,imod,ix,ipar,is)
              end do
            end do
          end do
        end do
      end if
    end do

    ! deallocate the field_id
    if (allocated(field_id)) deallocate(field_id)

    ! initialise the factors
    fnorm1d = 1.

    accumulated_normfactor = 1.
  else if (normalized) then

    if(non_linear) then
      call gkw_warn('Run is nonlinear - normalised=F is forced.')
      normalized = .false.
    end if
    ! initialise the normalisation factor
    fnorm1d = 1.
    accumulated_normfactor = 1.
  else

    ! initialise the normalisation factor.
    fnorm1d = 1.
    accumulated_normfactor = 1.

  end if
  
  
  if(non_linear .and. normalize_per_toroidal_mode) then
    nonlin_norm_fac = 1e-10
  else
    nonlin_norm_fac = 1.0
  endif

  ! done with the intialisation
  initialised = .true.

end subroutine normalise_init


!-------------------------------------------------------------------------------
!> Normalise the input fdis by a factor calculated from the fields.
!------------------------------------------------------------------------------
subroutine normalise_fdisi(fdis,nsolc)
  use control, only : normalized, normalize_per_toroidal_mode
  use dist, only : nregular_fields_end
  
  integer, intent(in) :: nsolc
  complex, intent(inout) :: fdis(nsolc)

  call calc_norm_factor(fdis, nregular_fields_end)
  ! perform the actual normalisation
  if (normalized .or. normalize_per_toroidal_mode) then
    call normalise_fdisk_with_fdisi(fdis,nsolc)
  end if

end subroutine normalise_fdisi


!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine calc_norm_factor(fdis,nregular_fields_end)
  use dist,         only : number_of_fields
  use mpicomms,     only : COMM_S_NE_X_NE
  use mpiinterface, only : mpiallreduce_sum_inplace
  use global,       only : r_tiny
  use grid,         only : parallel_s, parallel_x, nmod, nx, ns
  use geom,         only : ints
  use control,      only : normalized, normalize_per_toroidal_mode

  integer, intent(in) :: nregular_fields_end
  complex, intent(inout) :: fdis(nregular_fields_end)

  integer :: imod, j

  call calc_amp(fdis,nregular_fields_end)

  if (normalize_per_toroidal_mode) then
    ! Normalise with a separate normalisation factor for each toroidal
    ! mode.  Note that in each toroidal modes there may actually be
    ! several linear modes (see nmodes in mode module),
    ! if ikxspace > 1, because of the linear coupling due to the
    ! shear-periodicity at the parallel boundary conditions.
    
    ! calculate the normalisation factors
    amp_per_mode(:) = 0.
    ! sum the modulus square of the regular fields:
    do j = 1, number_of_fields*nx*ns
      do imod = 1, nmod
        amp_per_mode(imod) = amp_per_mode(imod) + abs(fdis(lindx(imod,j)))**2
      end do
    end do

    ! reduce sum over space and take sqrt()
    if (parallel_s .or. parallel_x) then
      ! sum over other processors responsible for other parts of fields
      call mpiallreduce_sum_inplace(amp_per_mode,nmod,COMM_S_NE_X_NE)
    end if

    ! make the norm factor  { |phi|^2 + |Apar|^2 + |Bpar|^2}
    amp_per_mode = amp_per_mode*ints(1)

    ! check if the value is too small;
    where (amp_per_mode < r_tiny)
      amp_per_mode = 1.
    elsewhere
      amp_per_mode = sqrt(amp_per_mode)
    end where

    fnorm1d = amp_per_mode
  else if (normalized) then
    ! Normalise with a single normalisation factor
    fnorm1d = amp
  else
    fnorm1d = 1.0
  end if

end subroutine calc_norm_factor

!----------------------------------------------------------------------------
!> Calculate this even if normalized == false
!> because it is needed to compute growth rates.
!----------------------------------------------------------------------------
subroutine calc_amp(fdis, nregular_fields_end)
  use dist, only : n_phi_start
  use grid, only : parallel_x, parallel_s
  use mpiinterface, only : mpiallreduce_sum_inplace
  use mpicomms, only : COMM_S_NE_X_NE
  use geom, only : ints
  use global, only : r_tiny
  integer, intent(in) :: nregular_fields_end
  complex, intent(inout) :: fdis(nregular_fields_end)

  ! first sum of local contribution
  amp = sum(abs(fdis(n_phi_start:nregular_fields_end))**2)

  ! reduce sum over processors
  if (parallel_s .or. parallel_x) then
    ! sum over other processors responsible for other parts of s
    call mpiallreduce_sum_inplace(amp,1,COMM_S_NE_X_NE)
  end if

  ! make the norm factor  { |phi|^2 + |Apar|^2 + |Bpar|^2}
  amp = amp*ints(1)

  ! check if the value is too small
  if (amp < r_tiny) then
    ! use 1.0 instead of small values
    amp = 1.
  else
    amp = sqrt(amp)
  endif
end subroutine calc_amp



!****************************************************************************
!> Normalise the input fdisk by a factor calculated previously, (possibly
!> using different data).
!> This routine is where the normalisation actually takes place.
!----------------------------------------------------------------------------
subroutine normalise_fdisk_with_fdisi(fdisk,nsolc)
  use control, only : normalized, normalize_per_toroidal_mode, non_linear
  use grid, only : nmod
  use dist, only : nf
  use global, only : r_huge

  integer, intent(in) :: nsolc
  !> as the fields are not needed any more, nsolc elements are enough
  complex, intent(inout) :: fdisk(nsolc)

  real :: one_over_factor
  real :: one_over_factor_per_mode(nmod)
  integer :: imod, j

  if (normalize_per_toroidal_mode) then
  
    if(non_linear) then
      one_over_factor_per_mode = nonlin_norm_fac / fnorm1d
    else
      one_over_factor_per_mode = 1. / fnorm1d
    endif
    
    do j = 1, nf/nmod
      do imod = 1, nmod
        fdisk(findx(imod,j)) = fdisk(findx(imod,j)) * one_over_factor_per_mode(imod)
      end do
    end do
    ! normalise the regular fields, the gyroavg fields, and the
    ! collisionop-related fields
    do j = 1, (nsolc-nf)/nmod
      do imod = 1, nmod
        fdisk(lindx(imod,j)) = fdisk(lindx(imod,j)) * one_over_factor_per_mode(imod)
      end do
    end do

    accumulated_normfactor = accumulated_normfactor * one_over_factor_per_mode

  else if (normalized) then
    ! perform on the whole array when a single value is used
    one_over_factor = 1. / fnorm1d(1)
    fdisk = fdisk * one_over_factor
    ! accumulate the normalisation factors, this can be handy for diagnostics.
    ! Eventually, accumulated_normfactor can become too tiny and
    ! looses its sense then!
    accumulated_normfactor = accumulated_normfactor * one_over_factor
  else
    ! do nothing if not normalised in any way
  end if

end subroutine normalise_fdisk_with_fdisi

!----------------------------------------------------------------------------
!>
!----------------------------------------------------------------------------
subroutine rotate_fdis(fdis)
  use grid, only : nmod, nx, nsp, ns, nmu, nvpar
  use control, only : nlphi, nlapar, nlbpar
  use dist, only : iphi, iapar, ibpar, ifdis
  use index_function, only : indx
  complex, intent(inout) :: fdis(:)
  complex :: rotate(nmod)
  integer :: imod, ix, is, i, k, j

  rotate = get_cmplx_rotation_factor()

  ! at the moment this rotates only the distribution and the
  ! not-gyroavged potentials.
  do imod=1, nmod
    do ix = 1, nx
      do is = 1, nsp
        do i = 1, ns
          do j = 1, nmu
            do k = 1, nvpar
              fdis(indx(ifdis,imod,ix,i,j,k,is)) = &
                 & fdis(indx(ifdis,imod,ix,i,j,k,is)) * rotate(imod)
            end do
          end do
        end do
      end do
    end do
  end do
  if(nlphi) then
    do imod=1, nmod
      do ix = 1, nx
        do is = 1, nsp
          do i = 1, ns
            do j = 1, nmu
              do k = 1, nvpar
                fdis(indx(iphi,imod,ix,i)) = fdis(indx(iphi,imod,ix,i)) * &
                   & rotate(imod)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
  if(nlapar) then
    do imod=1, nmod
      do ix = 1, nx
        do is = 1, nsp
          do i = 1, ns
            do j = 1, nmu
              do k = 1, nvpar
                fdis(indx(iapar,imod,ix,i)) = fdis(indx(iapar,imod,ix,i)) * &
                   & rotate(imod)
              end do
            end do
          end do
        end do
      end do
    end do
  end if
  if(nlbpar) then
    do imod=1, nmod
      do ix = 1, nx
        do is = 1, nsp
          do i = 1, ns
            do j = 1, nmu
              do k = 1, nvpar
                fdis(indx(ibpar,imod,ix,i)) = fdis(indx(ibpar,imod,ix,i)) * &
                   & rotate(imod)
              end do
            end do
          end do
        end do
      end do
    end do
  end if

end subroutine rotate_fdis

!----------------------------------------------------------------------------
!>
!----------------------------------------------------------------------------
function get_cmplx_rotation_factor()
  use control, only : non_linear, flux_tube
  use mpiinterface, only : mpiallreduce_maxloc, mpibcast, processor_number
  use mpiinterface, only : root_processor
  use dist, only : get_phi, phi, fdisi
  use mode, only : ixzero
  use global, only : r_tiny
  use grid, only : nmod, nx, ns

  complex, dimension(nmod) :: get_cmplx_rotation_factor

  integer :: imod
  integer :: ihelp
  integer, dimension(2) :: ihelp2
  real, dimension(2) :: phimax_local,phimax_global
  complex :: rotate

  get_cmplx_rotation_factor = 1.0

  if (.not.non_linear) then

    ! FIXME how many calls of get_phi are really necessary?
    call get_phi(fdisi,phi)

    do imod = 1, nmod

      ! Find the maximum of the potential and with (mpi) maxloc
      ! Find a (single) processor on which that point exists
      if (flux_tube) then
        ihelp=maxloc(abs(phi(imod,ixzero,1:ns)),1)
        ! beware using maxloc with array indices not starting at 1
        phimax_local(1)=abs(phi(imod,ixzero,ihelp))
        rotate=phi(imod,ixzero,ihelp)
      else
        ihelp2 = maxloc(abs(phi(imod,1:nx,1:ns)))
        !beware using maxloc with array indices not starting at 1
        rotate = phi(imod,ihelp2(1),ihelp2(2))
        phimax_local(1) = abs(rotate)
      end if

      ! use real(processor_number) for maxloc interface
      ! (could make a more Fortran friendly wrapper here)
      phimax_local(2)=real(1.*processor_number)
      call mpiallreduce_maxloc(phimax_local,phimax_global,1)
      ! convert back to integer for comparison
      ihelp=int(phimax_global(2))
      !Broadcast the point with maximum potential to all.
      call mpibcast(rotate,2,PROC=ihelp)
      !Avoid dividing by zero
      if (abs(rotate) < r_tiny) rotate = (1.,0.)

      ! norm it to 1
      rotate = rotate / abs(rotate)
      ! inverse the phase
      rotate = 1.0 / rotate

      get_cmplx_rotation_factor(imod) = rotate

    end do

    if (root_processor) then
      write(*,*) 
      write(*,*) 'mode structure: normalized (rotated in complex plane) by: '
      write(*,*) get_cmplx_rotation_factor
      write(*,*)
    end if
  end if


end function get_cmplx_rotation_factor


subroutine normalise_after_timestep()
  use dist, only : fdisi, fdisk, rhsk, nsolc
  use control, only : method, meth
  integer :: i

  call normalise_fdisi(fdisi,nsolc)

  if (method == 'EXP' .and. meth == 3) then
    ! normalize properly
    ! done after calculate_fields to have the new potential
    do i = 1, size(fdisk,2)
      call normalise_fdisk_with_fdisi(fdisk(:,i),nsolc)
    end do
    do i = 1, size(rhsk,2)
      call normalise_fdisk_with_fdisi(rhsk(:,i),nsolc)
    end do
  end if
end subroutine normalise_after_timestep


end module normalise
