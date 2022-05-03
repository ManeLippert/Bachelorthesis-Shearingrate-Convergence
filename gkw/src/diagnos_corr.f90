!------------------------------------------------------------------------------
!> The correlation function diagnostic
!------------------------------------------------------------------------------
module diagnos_corr

  implicit none

  private

  !>==== DESCRIPTION ==================================================
  !>
  !> The Diagnostic calulates the correlation function of the electric
  !> potential in radial and binomal direction.
  !>
  !> For spectral radius=true it uses the convolution theorem
  !> while it is calculated directly for spectral radius =false.
  !>
  !> Further a correlation length is calculated with 3 Definitions:
  !>  - The integral over the correlation function
  !>  - The 1/e point of the correlation function
  !>  - The first minimum of the correlation function (usefull in y direction)
  !>
  !>  For spectral radius=true one can also calculate the correlation length 
  !>  for each mode by seting imod_corr = (/imod ,xmod/)
  !>  it  will calculate the radial correlation length for the first ymod 
  !>  binormal modes and the binormal correlation for xmod radial modes
  !>  (note the radial modes start with -kmax and have the 0 mode at NX/2
  !>  the second half is mirrored to the first)
  !>  Correlation_length_x now has 3*ymod  and
  !>  Correlation_length_y 3*xmod values
  !>  the correlation functions are written imod/xmod lines for every timestep
  !>  
  !>   
  !>
  !>---- LIMITATIONS --------------------------------------------------
  !>
  !> imod_corr  is set to 1,1 for nonspectral. It should be possible to add 
  !> mode seperation here too, but is not yet implementated.
  !>
  !>===================================================================


  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output, final_output
  public :: output

  !> The general on/off switch for this diagnostic
  logical, save, public :: lcalc_corr
  !> (only for spectral_radius = T) Number of modes for which the data is
  !> calculated separately. The first number is the number of toroidal
  !> modes for which the correlation length is calculated, the second is
  !> the number of radial modes. If the values 1,1 are given, the
  !> correlation function will be averaged over all modes.
  integer,save, public :: imod_corr(2)
  !> The logical unit number, used to output data
  integer, save :: lun_corr_x=-1
  integer, save :: lun_corr_y=-1
  integer, save :: lun_corr_length_x=-1
  integer, save :: lun_corr_length_y=-1

  !> Variales and buffer for calculation
  real, save, allocatable :: buffer_x (:,:)
  real, save, allocatable :: buffer_y (:,:)
  complex, allocatable, save :: ax(:,:) , agx(:,:) , corr_compx(:) 

  !>power Spectrum
  complex, save, allocatable :: power_x(:)
  complex, save, allocatable :: power_y(:)
  real, save, allocatable :: power_real_x(:,:)
  real, save, allocatable :: power_real_y(:,:)

  !> Avaraged Phi
  real, save, allocatable :: avgphi(:)

  !> correlation funktion
  real, save, allocatable :: corr_fx(:,:)
  real, save, allocatable :: corr_fy(:,:)
  !>correlation length
  real, save, allocatable ::corr_l_x(:,:)
  real, save, allocatable ::corr_l_y(:,:)


contains

  subroutine set_default_nml_values()
    ! The diagnostic should be switched
    ! off by default.

    lcalc_corr = .false.
    imod_corr = (/ 1 ,1 /)
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_corr,1)
    call mpibcast(imod_corr,2)

  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_warn
    use grid,    only :  nmod ,n_x_grid
    use control, only : spectral_radius
    use mode,    only : mode_box
    if (.not.lcalc_corr) return

    if(imod_corr(1) > nmod) then
      call gkw_warn('diagnostic corr:: imod_corr(1) must be smaller than nmod')
      lcalc_corr = .false.
    end if

    if(imod_corr(2) > n_x_grid) then
      call gkw_warn('diagnostic corr:: imod_corr(2) must be smaller than n_x_grid')
      lcalc_corr = .false.
    end if

    if(.not.spectral_radius ) then
      if(imod_corr(1) /= 1 .and. imod_corr(2) /= 1 ) then
        imod_corr = 1
        call gkw_warn('imod_corr(:) is overwritten and set to 1')
      end if
    end if

    if(.not. mode_box) then
      call gkw_warn('diagnostic corr:: mode_box = true is needed for the correlation diagnostic')
      lcalc_corr = .false.
    end if

    if(nmod < 3) then
      call gkw_warn('diagnostic corr:: nmod > 2 is needed for a correlation length')
      lcalc_corr = .false.
    end if

    if(n_x_grid < 3) then
      call gkw_warn('diagnostic corr:: n_x_grid > 2 is needed for a correlation length')
      lcalc_corr = .false.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use diagnos_generic, only : attach_metadata_grid, LOCAL_DATA
    use grid, only : n_x_grid, nmod
    use control, only : spectral_radius
    use global, only : int2char, PHI_FIELD
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    logical, intent(inout) :: requirements(:,:)
    integer :: xsize
    integer :: ysize

    if (.not.lcalc_corr) return

    requirements(PHI_FIELD,LOCAL_DATA) = .true.

    if(root_processor) then
      xsize = n_x_grid
      if (spectral_radius) then
        xsize = fftsize(xsize)
      endif
      ysize = 2*nmod-2
      ysize = fftsize (ysize)     

      ! open logical units
      call open_real_lu('correlation_function_x', 'diagnostic/diagnos_corr', &
         & (/ xsize, imod_corr(1) /), ascii_fmt, lun_corr_x)
      call attach_metadata_grid(lun_corr_x, 'time', 'displacement', ascii_fmt) 
      call attach_metadata(lun_corr_x, description_key, &
         & 'X-correlation-function', ascii_fmt)
      if(imod_corr(1) == 1) then
        call attach_metadata(lun_corr_x, comments_key, &
         & 'radial correlation function integrated over s and y &
         & ', &
         & ascii_fmt)
      else
        call attach_metadata(lun_corr_x, comments_key, &
         & 'radial correlation function &
         & calculated separately for the '//int2char(imod_corr(1),2) &
         & //' binormal modes (according to the imod_corr(1) parameter) &
         & with the smallest wavenumbers.', &
         & ascii_fmt)
      end if
      call attach_metadata(lun_corr_x, phys_unit_key, not_avail, ascii_fmt)

      call open_real_lu('correlation_function_y', 'diagnostic/diagnos_corr', &
         & (/ ysize, imod_corr(2) /), ascii_fmt, lun_corr_y)
      call attach_metadata_grid(lun_corr_y, 'time', 'Displacement', ascii_fmt) 
      call attach_metadata(lun_corr_y, description_key, &
         & 'Y-correlation-function', ascii_fmt)
      if(imod_corr(1) == 1) then
        call attach_metadata(lun_corr_y, comments_key, &
         & 'binormal correlation function integrated over s and x &
         & ', &
         & ascii_fmt)
      else
        call attach_metadata(lun_corr_y, comments_key, &
         & 'binormal correlation function &
         & calculated separately for the '//int2char(imod_corr(2),2) &
         & //' radial modes (according to the imod_corr(2) parameter) &
         & with the smallest wavenumbers.', &
         & ascii_fmt)
      end if
      call attach_metadata(lun_corr_y, phys_unit_key, not_avail, ascii_fmt)


      if(imod_corr(1) > 1) then
        call open_real_lu('correlation_length_x', 'diagnostic/diagnos_corr', &
           & (/ 3,imod_corr(1) /), ascii_fmt, lun_corr_length_x)
        call attach_metadata_grid(lun_corr_length_x, &
           & 'time', '3 correlation lengths', 'krho(1:imod_corr(1))', ascii_fmt)
      else
        call open_real_lu('correlation_length_x', 'diagnostic/diagnos_corr', &
         & (/ 3 /), ascii_fmt, lun_corr_length_x)
        call attach_metadata_grid(lun_corr_length_x, &
           & 'time', '3 correlation lengths', ascii_fmt)
      end if
      call attach_metadata(lun_corr_length_x, description_key, &
         & 'correlation length in X  direction', ascii_fmt) 
      call attach_metadata(lun_corr_length_x, comments_key, &
         & 'The three definitions of the correlation length are &
         & (i) integral over the correlation funktion &
         & (ii) 1/e point of the correlation function &
         & (iii) first Minimum of the correlation function. &
         & ', &
         & ascii_fmt)
      call attach_metadata(lun_corr_length_x, phys_unit_key, not_avail, ascii_fmt)
      
      if(imod_corr(2) > 1) then
        call open_real_lu('correlation_length_y', 'diagnostic/diagnos_corr', &
           & (/ 3,imod_corr(2) /), ascii_fmt, lun_corr_length_y)
        call attach_metadata_grid(lun_corr_length_y, &
           & 'time', '3 correlation lengths', 'kxrh(1:imod_corr(1))', ascii_fmt)
      else
        call open_real_lu('correlation_length_y', 'diagnostic/diagnos_corr', &
           & (/ 3 /), ascii_fmt, lun_corr_length_y)
        call attach_metadata_grid(lun_corr_length_y, &
           & 'time', '3 correlation lengths', ascii_fmt)
      end if
      call attach_metadata(lun_corr_length_y, description_key, &
         & 'correlation length in Y direction', ascii_fmt) 
      call attach_metadata(lun_corr_length_y, comments_key, &
         & 'The three definitions of the correlation length are &
         & (i) integral over the correlation funktion &
         & (ii) 1/e point of the correlation function &
         & (iii) first Minimum of the correlation function. &
         & ', &
         & ascii_fmt)
      call attach_metadata(lun_corr_length_y, phys_unit_key, not_avail, ascii_fmt)

    end if

  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
    use control,      only : spectral_radius
    use grid,    only : nmod, n_x_grid, nx
    use general, only : gkw_abort 
    use control,      only : spectral_radius
    integer :: ierr
    integer :: xsize
    integer :: ysize

    if (.not.lcalc_corr) return
    !allocate correlation function length = the next power of 2 after nx, ny
    xsize = n_x_grid
    if (spectral_radius) then
      xsize= fftsize (xsize)
    endif
    ysize = 2*nmod-2
    ysize = fftsize (ysize)
    

    !allocate(corr_fx(xsize),stat=ierr)
    allocate(corr_fx(xsize, imod_corr(1)),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(corr_fy(ysize, imod_corr(2)),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')  

    allocate(corr_l_x(3,imod_corr(1)),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate') 
    allocate(corr_l_y(3,imod_corr(2)),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate') 

    allocate(buffer_x(n_x_grid,nmod), stat = ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate buffer_x')
    allocate(buffer_y(nmod,n_x_grid), stat = ierr)  
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate buffer_y')


    allocate(power_real_x(n_x_grid,nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')
    allocate(power_real_y(nmod,n_x_grid),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')

    allocate(avgphi(n_x_grid),stat=ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')


    allocate(power_x((xsize/2+1)), stat = ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')       
    allocate(power_y((ysize/2+1)), stat = ierr)
    if (ierr /= 0) call gkw_abort('diagnostic :: could not allocate')

    if (.not. spectral_radius) then
      allocate(ax(nx,nmod), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate ax')

      allocate(agx(n_x_grid,nmod), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate agx')

      allocate(corr_compx(n_x_grid), stat = ierr)
      if (ierr /= 0) call gkw_abort('could not allocate filxy')
    end if

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu , ascii_fmt
    use mpiinterface, only : root_processor
    ! deallocate all arrays of this diagnostic
    !if(allocated(local_array_name)) deallocate(local_array_name)

    if (.not.lcalc_corr) return
    if(.not. root_processor) return
    ! be nice and close all logical units
    call close_lu(lun_corr_x, ascii_fmt)
    call close_lu(lun_corr_y, ascii_fmt)
    call close_lu(lun_corr_length_x, ascii_fmt)
    call close_lu(lun_corr_length_y, ascii_fmt)


    if(allocated(corr_fx)) deallocate(corr_fx)
    if(allocated(corr_fy)) deallocate(corr_fy)
    if(allocated(buffer_x)) deallocate(buffer_x)
    if(allocated(buffer_y)) deallocate(buffer_y)    
    if(allocated(corr_l_x)) deallocate(corr_l_x)
    if(allocated(corr_l_y)) deallocate(corr_l_y)
    if(allocated(power_x)) deallocate(power_x)
    if(allocated(power_y)) deallocate(power_y)
    if(allocated(power_real_x)) deallocate(power_real_x)
    if(allocated(power_real_y)) deallocate(power_real_y)
    if(allocated(ax)) deallocate(ax)
    if(allocated(agx)) deallocate(agx)
    if(allocated(corr_compx)) deallocate(corr_compx)
    if(allocated(avgphi)) deallocate(avgphi)
  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

  end subroutine initial_output


  subroutine final_output(number)
    integer, intent(in), optional :: number

    if (.false.) write(*,*) number ! to keep the compiler quiet.

  end subroutine final_output

  !--------------------------------------------------------------------
  !> the function gives the size for a FFT array  for zero padding 
  !--------------------------------------------------------------------
  function fftsize(grid) result(outgrid)
    integer, intent(inout) :: grid
    integer :: outgrid
    real :: xsize
    integer :: i

    ! look for an integer which is at least larger than grid by a factor,
    ! and a an efficient grid size for the FFT
    ! (at the moment a power of two)
    i = 1
    xsize = grid*1.3
    do while (xsize > 2)
      xsize = xsize / 2
      i = i+1
    end do
    outgrid = 2**i

  end function fftsize

  
  !--------------------------------------------------------------------
  !> largestep calculates the correlation function at every timestep
  !> and calls the routine to calculate the correlation length
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use grid,         only : nmod, nx, ns, n_x_grid, ls, n_s_grid
    use grid,         only : proc_subset
    use dist,         only : phi
    use control,      only : spectral_radius
    use fft,          only : four1D_real, FFT_INVERSE
    use geom,         only : ints
    use mpiinterface, only : mpiallreduce_sum, gather_array, root_processor, mpiallreduce_sum_inplace
    use mpicomms,     only : COMM_S_NE,  COMM_DUMMY , COMM_X_NE , COMM_S_NE_X_NE
    use mode,         only : lyn, lxn, ixzero
    use diagnos_generic, only : parseval_correction
    !> counting indices
    integer ::  imod, ix, is ,xmod
    integer :: r
    integer :: fftsize_x
    integer :: fftsize_y
    complex :: fac

    if (.not.lcalc_corr) return

    if (spectral_radius) then
      !> X-correlation Function
      !----------------------------

      ! Calculate the x power spectrum   
      buffer_x = 0.
      do ix = 1, n_x_grid 
        do is = 1, n_s_grid
          if(proc_subset(0,is,0,0,0)) then !For parallel_s
            do imod = 1, nmod    
              buffer_x(ix,imod) = buffer_x(ix,imod) + abs(phi(imod,ix,ls(is)))**2 *ints(is) * parseval_correction(imod)
            end do
          end if
        end do
      end do

      ! finish of the fluxsurface average by a mpireduce over s direction

      power_real_x = 0.
      call mpiallreduce_sum(buffer_x,power_real_x,n_x_grid,nmod,COMM_S_NE)


      do imod=1,imod_corr(1)
        !for standart case sum over the modes otherwise the modes are treated seperatly 
        if (imod_corr(1) .eq. 1) then
          power_real_x(:,1) = 0.
          do ix = 1, n_x_grid
            power_real_x(ix,1) = sum(power_real_x(ix,:),1)
          end do
        end if
        ! power_real_x now contains the x power-spectrum of phi
        ! Note that the coefficient belonging to the 0 wavenumber is in the middle of the array
        ! (see also ixzero) 

        if(root_processor) then

          !> pad the FFT array with zeros
          power_x(:) = 0. 
          fftsize_x = n_x_grid
          fftsize_x = fftsize(fftsize_x)
          

          !> inverse Fouriertransformation for the correlation funktion
          corr_fx(:,imod) = 0.

          power_x(1:ixzero) = power_real_x(ixzero:n_x_grid,imod)
          
          call four1D_real(corr_fx(:,imod),power_x,FFT_INVERSE)

          !> calculate Correlation length

          call correlation_length(corr_fx(:, imod),fftsize_x,lxn,corr_l_x(:,imod))
        end if
      end do


      !> Y-correlation Function
      !----------------------------

      buffer_y = 0.
      !>integrate |phi(ky)|^2 over s for the power spectrum

      do imod = 1, nmod
        do is = 1, n_s_grid 
          if(proc_subset(0,is,0,0,0)) then !For parallel_s
            do ix = 1, n_x_grid             
              buffer_y(imod,ix) = buffer_y(imod,ix) + abs(phi(imod,ix,ls(is)) )**2 *ints(ls(is))
            end do
          end if
        end do
      end do
      power_real_y = 0.
      call mpiallreduce_sum(buffer_y,power_real_y,nmod,n_x_grid,COMM_S_NE)

      !calculate the binormal average of phi which equals the 0 mode

      avgphi = 0.
      do is = 1, n_s_grid
        if(proc_subset(0,is,0,0,0)) then !For parallel_s
          do ix = 1, n_x_grid
            avgphi(ix) = avgphi(ix) + abs(phi(1,ix,ls(is)) )**2 *ints(ls(is))
          end do
        end if
      end do
      !> finish of the fluxsurface average by a mpireduce over s direction

      call mpiallreduce_sum_inplace(avgphi,n_x_grid,COMM_S_NE)


      do xmod=1,imod_corr(2) 

        if ( root_processor) then
          !for standart case sum over x otherwise the xmodes are treated seperatly 
          if (imod_corr(2) .eq. 1) then
            do ix = 2, n_x_grid  
              do imod = 1, nmod    
                power_real_y(imod,1) =  power_real_y(imod,1) + power_real_y(imod,ix) 
              end do
              avgphi(1) =avgphi(1) +avgphi(ix)
            end do
          end if


          !> zero padding the FFT array
          power_y(:) = 0. 
          fftsize_y = 2*nmod-2
          fftsize_y= fftsize(fftsize_y)

          !> inverse Fouriertransformation for the correlation funktion
          corr_fy(:,xmod) = 0.
          power_y(1:nmod) = power_real_y(1:nmod,xmod)
          call four1D_real(corr_fy(:,xmod),power_y, FFT_INVERSE)

          do imod=1,fftsize_y     
            corr_fy(imod,xmod) = corr_fy(imod,xmod) - avgphi(xmod)
          end do

          !>  corelation length  
          call correlation_length(corr_fy(:,xmod),fftsize_y,lyn,corr_l_y(:,xmod))
          
        end if
      end do

    else
      !------------------------------------------------------------------------
      ! nonspectral case
      !------------------------------------------------------------------------
      !> X-correlation Function
      !----------------------------
      agx = 0.
      ax = 0.
      buffer_x = 0.    
      corr_compx = 0.
      do is = 1, ns        
        do imod = 1, nmod
          do ix = 1, nx
            ax(ix, imod) = phi(imod,ix,is)
          end do
        end do

        ! gather in x direction
        call gather_array(agx, n_x_grid, nmod, &
          & ax, nx, nmod, &
          & COMM_X_NE, COMM_DUMMY,.true.)

        if(proc_subset(1,0,0,0,0)) then
          do imod = 2, nmod
            do r = 1 , n_x_grid
              do ix = 1, n_x_grid
                fac = ints(is) * parseval_correction(imod)
                if ( (r + ix) < (n_x_grid +2)) then
                  corr_compx(r) = corr_compx(r)+CONJG(agx(ix,imod))*(agx(ix+r-1,imod)) * fac
                else
                  corr_compx(r) = corr_compx(r)+CONJG(agx(ix,imod))* (agx(ix+r-n_x_grid,imod))* fac
                end if
              end do
            end do
          end do
        end if

      end do

      call mpiallreduce_sum_inplace(corr_compx,n_x_grid,COMM_S_NE)


      if(root_processor) then
        corr_fx = 0.
        corr_fx(1:n_x_grid,1) = real(corr_compx)

        !>  corelation length
        call correlation_length(corr_fx(:,1),n_x_grid,lxn,corr_l_x(:,imod_corr(1)))

      end if


      !> Y-correlation Function
      !----------------------------
      
      buffer_y = 0.

      do is = 1, ns
        do imod = 1, nmod         
          do ix = 1, nx 
            buffer_y(imod,1) = buffer_y(imod,1) + abs(phi(imod,ix,is))**2 *ints(is) !* parseval_correction(imod)
          end do
        end do
      end do

      !> finish of the fluxsurface average by a mpireduce over s direction
      power_real_y = 0.
      call mpiallreduce_sum(buffer_y(1:nmod,1),power_real_y(1:nmod,1) ,nmod,COMM_S_NE_X_NE)

      ! do xmod=1,imod_corr(2)
      xmod=1
      avgphi = 0.
      do is = 1, n_s_grid
        if(proc_subset(0,is,0,0,0)) then !For parallel_s
          do ix = 1, nx
            avgphi(1) = avgphi(1) + abs(phi(1,ix,ls(is)) )**2 *ints(ls(is))
          end do
        end if
      end do
      !> finish of the fluxsurface average by a mpireduce over s direction

      call mpiallreduce_sum_inplace(avgphi,n_x_grid,COMM_S_NE_X_NE)


      if(root_processor) then

        !> zero padding the FFT array
        power_y(:) = 0. 
        fftsize_y = 2*nmod-2
        fftsize_y = fftsize(fftsize_y)
        

        !> inverse Fouriertransformation for the correlation funktion
        corr_fy = 0.
        power_y(1:nmod) = power_real_y(1:nmod,1)
        call four1D_real(corr_fy(:,xmod),power_y,FFT_INVERSE)

        do imod=1,fftsize_y     
          corr_fy(imod,xmod) = corr_fy(imod,xmod) - avgphi(xmod)
        end do

        !>  corelation length
        
        call correlation_length(corr_fy(:,xmod),fftsize_y,lyn,corr_l_y(:,xmod))
        
      end if

      ! end do !xmod
      
    end if !end nonspectral

  end subroutine calc_largestep

  !--------------------------------------------------------------------
  !> The routine calculates the correlation length from a given correlation function
  !> Input : the correlation function (corr), number of points of the correlation function (length)
  !> box size (dx), array to input the correlation length (corrlength) atm 3 values.
  !> Also does the normalisation of the correlation function.
  !-------------------------------------------------------------------
  subroutine correlation_length(corr, length,dx,corrlength)
    real, dimension(:), intent(inout) :: corr 
    integer, intent(in) :: length 
    real, intent(in) :: dx 
    real, dimension(:), intent(inout) :: corrlength
    integer :: ix
    real :: m_e

    m_e = exp(1.0)

    !>Normalisation of correlation_funktion
    if (corr(1) < tiny( corr(1))) then
      !set corrlength = half boxsize if the reduced potential is zero
      corrlength(1)= dx/2
      corrlength(2)= dx/2
      corrlength(3)= dx/2
    else
      corr =corr/corr(1)

      !> Integration over correlation function
      corrlength(1)=0
      do ix = 1, int((length+1)/2)
        corrlength(1) = corrlength(1)+corr(ix)*dx/length
      end do

      !> 1/e point
      ix = 1
      do while ( corr(ix) > 1/m_e  .and. ix < length )
        ix =ix+1
      end do
      ! linear interpolation of the 1/e point
      if (ix < (length+1)/2) then
        corrlength(2) =( ((1/m_e-corr(ix-1))/(corr(ix)-corr(ix-1)))  + ix - 2 ) *dx/length
      !set corrlength = half boxsize if there is no 1/e point
      else
        corrlength(2) = dx/2
      end if

      !> Minimum point
      ix = 2
      do while ( ix < length-2 .and. corr(ix) > corr(ix+1)   )
        ix =ix+1
      end do
      !interpolarion of the minimum
      if (ix < (length+1)/2) then
        corrlength(3) = (ix-(corr(ix+1)-corr(ix-1))/(2*corr(ix+1)+2*corr(ix-1)-4*corr(ix)) )*dx/length
        !set corrlength = half boxsize if there is no minimum
      else
        corrlength(3) = dx/2
      end if

      !> Integration over absolutum of correlation function
      !corrlength(2)=0
      !do ix = 1, int((length+1)/2)
      ! corrlength(4) = corrlength(4)+abs(corr(ix))*dx/length
      !end do

      !> Integral over r^2 abs
      !corrlength(3)=0
      !do ix = 1, int((length+1)/2)
      !  corrlength(5) = corrls
      !length(5)+ix*ix*abs(corr(ix))*dx/length
      !end do
      ! corrlength(5) = sqrt(corrlength(5) /corrlength(4)) *dx/length
    end if
  end subroutine correlation_length

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()
    use mpiinterface, only : root_processor
    use io, only : append_chunk, xy_fmt, ascii_fmt

    if (.not.lcalc_corr) return

    call calc_largestep

    if(root_processor) then
      
      call append_chunk(lun_corr_x, corr_fx, xy_fmt, ascii_fmt)  
      call append_chunk(lun_corr_y, corr_fy, xy_fmt, ascii_fmt)
      if(imod_corr(1) > 1) then
        call append_chunk(lun_corr_length_x, corr_l_x, xy_fmt, ascii_fmt)
      else
        call append_chunk(lun_corr_length_x, corr_l_x(:,1), xy_fmt, ascii_fmt)
      end if
      if(imod_corr(2) > 1) then
        call append_chunk(lun_corr_length_y, corr_l_y, xy_fmt, ascii_fmt)
      else
        call append_chunk(lun_corr_length_y, corr_l_y(:,1), xy_fmt, ascii_fmt)
      end if
      
    end if
  end subroutine output

end module diagnos_corr
