!------------------------------------------------------------------------------
!> Energy transfer diagnostic. Perhaps temporary 
!------------------------------------------------------------------------------
module diagnos_eng

  implicit none 

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output


  private 

  !> switch for a simple diagnostic on the mode energy (in the potential only)
  logical, save, public :: lmode_energy 


  !> Temporary storage for the distribution function and its time derivative 
  complex, save, allocatable :: fstore(:), dfdt(:) 

  !> zfs contains the time derivative of abs(phi)^2
  real,    save, allocatable :: zfs(:), duma(:) 

  ! the file identifier 
  integer, save :: lun1, lun2  

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    lmode_energy = .false.

    
  end subroutine set_default_nml_values


!------------------------------------------------------------------------------
!> Change in mode energy 
!------------------------------------------------------------------------------
subroutine mode_energy

  use dist,             only : fdisi, nsolc, iphi, ifdis
  use grid,             only : nmod, nx, ns, parallel_s, parallel_x, nmu 
  use grid,             only : nvpar, nsp 
  use index_function,   only : indx 
  use io,               only : xy_fmt, ascii_fmt, append_chunk
  use non_linear_terms, only : add_non_linear_terms
  use mpicomms,         only : COMM_S_NE_X_NE
  use mpiinterface,     only : mpiallreduce_sum
  use velocitygrid,     only : intvp, intmu 
  use geom,             only : bn 

  integer imm, imn, i, imod, ix, j, k, is 
  complex :: dndt 

  ! outer loop over the toroidal mode number. This mode number is what is 
  ! kept in the ExB velocity and f 
  do imm = 1, nmod; do imn = 1, nmod 

    ! copy the distribution function 
    do i = 1, nsolc 
      fstore(i) = fdisi(i) 
      dfdt(i)   = (0.,0.)
    end do 

    ! switch off all but the imm mode 
    do imod = 1, nmod; do ix = 1, nx; do i = 1, ns 
      if (imod .ne. imm) fstore(indx(iphi,imod,ix,i)) = (0.,0.) 
      if (imod .ne. imn) then 
        do j = 1, nmu; do k = 1, nvpar; do is = 1, nsp 
          fstore(indx(ifdis,imod,ix,i,j,k,is)) = (0.,0.)
        end do; end do; end do; 
      endif 
    end do; end do; end do 

    ! THIS IS IMPRACTICALLY SLOW to call inside this loop
    ! call the nonlinear time integrator 
    call add_non_linear_terms(fstore,dfdt)

    ! solve for the fields 
    !if (spectral_radius) then 
    !  call calculate_fields(dfdt)
    !else 
    !  call field_solve_nonspec_wrap(dfdt,0.)
    !endif 

    ! calculate the change in amplitude of the field 
    do imod = 1, nmod 
      zfs(imod) = 0. 
      duma(imod) = 0. 
      do ix = 1, nx; do i = 1, ns 

        ! calculate the density (conjugate) and the time derivative 
        dndt = (0.,0.) 
        do j = 1, nmu; do k = 1, nvpar 
          dndt = dndt + intvp(i,j,k,1)*intmu(j)*bn(ix,i)   *   & 
               & conjg(fdisi(indx(ifdis,imod,ix,i,j,k,1))) *   &
               & dfdt(indx(ifdis,imod,ix,i,j,k,1))    
        end do; end do 
        zfs(imod) = zfs(imod) + real(dndt)
      end do; end do 
    end do 

    ! communication 
    if (parallel_s .or. parallel_x) then
      ! sum over other processors responsible for other parts of feilds
      call mpiallreduce_sum(zfs,duma,nmod,COMM_S_NE_X_NE)
      do imod = 1, nmod 
        zfs(imod) = duma(imod) 
      end do
    end if

    ! write the result on file 
    call append_chunk(lun1, zfs, xy_fmt, ascii_fmt)

  end do; end do  ! loop over toroidal mode of the ExB velocity and f 

  do imod = 1, nmod 
    zfs(imod) = 0. 
    do ix = 1, nx; do i = 1, ns 

      ! calculate the density (conjugate) and the time derivative 
      dndt = (0.,0.) 
      do j = 1, nmu; do k = 1, nvpar 
        dndt = dndt + 0.5 * intvp(i,j,k,1) * intmu(j) * bn(ix,i) *     & 
             & abs(fdisi(indx(ifdis,imod,ix,i,j,k,1)))**2
      end do; end do 
      zfs(imod) = zfs(imod) + real(dndt)

    end do; end do 
  end do 
  call append_chunk(lun2, zfs, xy_fmt, ascii_fmt)

end subroutine mode_energy


!------------------------------------------------------------------------------
!> close the files, deallocate the arrays 
!------------------------------------------------------------------------------
subroutine finalize()
  use io, only : close_lu, ascii_fmt
  use mpiinterface, only : root_processor

  if (.not.lmode_energy) return
  
  if (root_processor) then
    call close_lu(lun1, ascii_fmt)
    call close_lu(lun2, ascii_fmt)

    ! deallocate the help arrays 
    if (allocated(fstore)) deallocate(fstore)
    if (allocated(dfdt))   deallocate(dfdt) 
    if (allocated(zfs))    deallocate(zfs)
    if (allocated(duma))   deallocate(duma)
  end if

end subroutine finalize

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lmode_energy, 1) 

  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use grid, only : nmod
    use mpiinterface, only : root_processor
    use io, only : open_real_lu, ascii_fmt
    use global, only : dotdat
    use control, only : io_legacy
    logical, intent(inout) :: requirements(:,:)

    ! To keep the compiler quiet, as the array is not used so far (in this routine).
    if (.false.) write(*,*) requirements

    if (.not. lmode_energy) return
    if(root_processor) then
      ! Open the file for output 
      call open_real_lu(dotdat('mode_energy',io_legacy), 'diagnostic/diagnos_eng', (/ nmod /), &
         & ascii_fmt, lun1)
      call open_real_lu(dotdat('entro',io_legacy), 'diagnostic/diagnos_eng',(/ nmod /), &
         & ascii_fmt, lun2)
    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use dist, only : nsolc
    use grid, only : nmod
    use general, only : gkw_abort
    integer :: ierr

    if (.not. lmode_energy) return
    
    ! allocate the arrays 
    ierr = 0  
    allocate(fstore(nsolc), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate fstore')
    allocate(dfdt(nsolc), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate dfdt')
    allocate(zfs(nmod), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate zfs') 
    allocate(duma(nmod), stat = ierr) 
    if (ierr /= 0) call gkw_abort('Could not allocate zfs') 

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine initial_output()

  end subroutine initial_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine final_output()

  end subroutine final_output

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output()
    
    ! The energy in the modes
    if (lmode_energy) call mode_energy
  end subroutine output

  
end module diagnos_eng
