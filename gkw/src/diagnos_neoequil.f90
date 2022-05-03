!------------------------------------------------------------------------------
!> Diagnostics related to the reader of NEO equilibria
!> 
!> DESCRIPTION
!> 
!> Module containing the diagnostics on the background, neoclassical
!> distribution function as calculated from NEO.  Calculated are:
!>
!> Neoclassical flows, bootstrap current and perturbed density and 
!> potential profiles
!>
!> Also output is the reconstructed distribution function in GKW velocity
!> space
!>
!>
!> LIMITATIONS: 
!>
!> Only currently useful for Neoclassical background
!>
!------------------------------------------------------------------------------
module diagnos_neoequil

  implicit none

  private
  
  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output
  public :: output, screen_output

  !> The general on/off switch for this diagnostic
  logical, save, public :: lcalc_neoequil = .true.

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    use linear_terms, only : lneo_equil, lneo_equil_switch_default

    lcalc_neoequil = ( lneo_equil .or. lneo_equil_switch_default )

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_neoequil,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general, only : gkw_warn
    use neoequil, only : neo_nsp
    use grid, only : nsp

    if(.not. lcalc_neoequil) return

    if(neo_nsp /= nsp) then
      call gkw_warn('diagnos_neoquil can not really deal with trace species &
         & so far')
      !lcalc_neoequil = .false.
    end if

  end subroutine check

  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init()
    use mpiinterface, only : root_processor

    if(.not. lcalc_neoequil) return
    
    if(root_processor) then
      ! Open logical units here.
      ! A PROPOSAL FOR NAMING RULES OF NEW LOGICAL UNITS: Note that the logical
      ! unit name does not
      ! contain the .dat, .bin or any other suffix.  Moreover the
      ! logical unit name should not contain dots (because in
      ! scripting languages the data is likely to be stored in a
      ! structure and there the dot is an operator and needs to be converted).
    end if
    
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()

    if(.not. lcalc_neoequil) return
    
  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()

    if(.not. lcalc_neoequil) return

  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()
    use neoequil, only : neo_n_species_incl_dup
    use grid, only : number_of_species
    use general, only : gkw_warn

    if(.not. lcalc_neoequil) return

    if(neo_n_species_incl_dup == number_of_species) then
      call gkw_warn('diagnos_neoquil can only deal with the simple case so far')
      call neof_diagnose()
      call output_reconstr_distr_pot()
    end if
    call output_gradients()

  end subroutine initial_output

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()

  end subroutine output

  !------------------------------------------------------------------------------
  !> 
  !------------------------------------------------------------------------------
  subroutine neof_diagnose

    use grid,      only : nvpar, ns, nmu, n_s_grid, nsp
    use grid,      only : nmu, number_of_species
    use geom,      only : bn, rfun, signB, bt_frac, ints 
    use velocitygrid,   only : vpgr, mugr, intvp, intmu
    use io,     only : get_free_file_unit
    use mpiinterface, only : root_processor
    use mpiinterface, only : mpiallreduce_sum_inplace, gather_array
    use mpiinterface, only : mpireduce_sum_inplace
    use components,   only : vthrat, mas, de, fp
    use dist,      only : fmaxwl
    use neoequil, only : neof, mid_r_point, gradneof, rrefoa, neophi
    use neoequil, only : neo_eps, neo_nsp
    use mpicomms,     only : COMM_S_NE,COMM_VPAR_NE_MU_NE
    use mpicomms,     only : COMM_SP_NE, COMM_S_NE
    use mpicomms,     only : COMM_SP_NE_S_NE, COMM_X_EQ
    use linear_terms, only : dmaxwel

    use grid, only : proc_subset
    use io, only : output_array, xy_fmt, ascii_fmt
    use general, only : gkw_abort
    
    real :: dum2
    real :: commf,dum3, commfgrad, commfgrad2, dum2grad
    real :: fluxsum2,rmean,fluxsum3,fluxsum2grad
    real :: velele
    real :: phisum(5),gradphi
    integer :: i,j,k,is,ix,file_unit,ixneo,iix
    real :: den_local(neo_nsp, ns), flow_local(neo_nsp, ns)
    real, allocatable :: den_global(:, :), flow_global(:, :), fluxsum1(:)
    real :: densum_local(ns),flowsum_local(ns)
    integer :: ierr

    ! pick the right subset to save memory
    if(proc_subset(0,0,0,0,0)) then
      allocate(den_global(number_of_species, n_s_grid), stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_neoequil :: could not allocate den_global')
      allocate(flow_global(number_of_species, n_s_grid),  stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_neoequil :: could not allocate flow_global')
      allocate(fluxsum1(number_of_species),  stat=ierr)
      if (ierr /= 0) &
         & call gkw_abort('diagnos_neoequil :: could not allocate fluxsum1')
    else
      allocate(den_global(1,1))
      allocate(flow_global(1,1))
      allocate(fluxsum1(1))
    end if
    
    if(root_processor)then
      call get_free_file_unit(file_unit)
      open(file_unit, status='unknown',action='write', file='neof_flows.dat')  
    endif

    if(root_processor)write(file_unit,*) 'All neoequil outputs are in GKW &
       & normalised units'
    
    !Only the middle neo flux tube is output.
    ixneo = mid_r_point
    ix = 1
    fluxsum1 = 0.E0
    fluxsum2 = 0.E0
    fluxsum2grad=0.E0
    fluxsum3 = 0.E0
    rmean = 0.E0
    den_local = 0.E0
    flow_local = 0.E0
    phisum = 0.E0
    do i = 1,ns 
      do iix = 1,5
        phisum(iix) = phisum(iix) + neophi(iix,i)**2*ints(i)
      end do
      neo_species: do is=1,neo_nsp
        dum2 = 0.E0
        dum3 = 0.E0
        dum2grad=0.E0
        do j = 1,nmu
          do k = 1,nvpar
            commf = fmaxwl(ix,i,j,k,is)*neof(ixneo,is,i,j,k)
            commfgrad = fmaxwl(ix,i,j,k,is)*gradneof(is,i,j,k)
            !Removes the density gradient part as it is in the denominator of the expression
            commfgrad2 =  (dmaxwel(ix,i,j,k,is,0)-fp(ix,is))*fmaxwl(ix,i,j,k,is)* &
               & neof(ixneo,is,i,j,k)

            velele = intvp(i,j,k,is)*intmu(j)*bn(ix,i)

            !Local density per species         
            den_local(is,i) = den_local(is,i) + commf*velele

            !Local parallel flow per species
            flow_local(is,i) = flow_local(is,i) + vthrat(is)*vpgr(i,j,k,is)*&
               & commf*velele

            !Toroidal diamagnetic flow  
            dum2 = dum2 + commf*vthrat(is)*vpgr(i,j,k,is)*Rfun(ix,i)*&
               & bt_frac(ix,i)*signB*velele
            dum2grad = dum2grad +  commfgrad*vthrat(is)*vpgr(i,j,k,is)* &
               & Rfun(ix,i)*bt_frac(ix,i)*signB*velele
            dum2grad = dum2grad -  commfgrad2*vthrat(is)*vpgr(i,j,k,is)* &
               & Rfun(ix,i)*bt_frac(ix,i)*signB*velele
            dum3 = dum3 + commf*vthrat(is)*sqrt(2*mugr(j)*bn(ix,i))* &
               & sqrt(1-bt_frac(ix,i)**2)*velele/Rfun(ix,i)
          enddo
        enddo

        fluxsum1(is) = fluxsum1(is) + bn(ix,i)*flow_local(is,i)*ints(i)
        fluxsum2 = fluxsum2 + mas(is)*dum2*ints(i)
        fluxsum2grad = fluxsum2grad + mas(is)*dum2grad*ints(i)
        fluxsum3 = fluxsum3 + dum3*ints(i)
        rmean = rmean + mas(is)*de(ix,is)*ints(i)*Rfun(ix,i)**2
      enddo neo_species

    enddo
    
    if(neo_nsp > 0) then
      call mpireduce_sum_inplace(den_local, shape(den_local), COMM_VPAR_NE_MU_NE)
      call mpireduce_sum_inplace(flow_local, shape(flow_local), COMM_VPAR_NE_MU_NE)
    end if
   
    ! FIXME these gathers -> gatherv or append_chunk like for
    ! gradients, for neo_nsp /= nsp
    call gather_array(den_global,number_of_species, n_s_grid, &
       & den_local,nsp,ns, &
       & COMM_SP_NE, COMM_S_NE)
    write (*,*) "gather 2nd"
    call gather_array(flow_global,number_of_species, n_s_grid, &
       & flow_local,nsp,ns, &
       & COMM_SP_NE, COMM_S_NE)

    if(root_processor) then
      call output_array('neof_den_profiles', 'diagnostic/diagnos_neoequil', &
         & den_global, 'F', &
         & xy_fmt, ascii_fmt)
      call output_array('neof_flow_profiles', 'diagnostic/diagnos_neoequil', &
         & flow_global, 'F', &
         & xy_fmt, ascii_fmt)
    end if

    if(proc_subset(0,0,1,1,0)) then
      densum_local = 0.0
      flowsum_local = 0.0
      do is = 1, neo_nsp
        do i = 1, ns
          densum_local(is) = densum_local(is) + den_local(is,i)*ints(i)
          flowsum_local(is) = flowsum_local(is) + flow_local(is,i)*ints(i)
        end do
      end do
      call mpireduce_sum_inplace(densum_local,shape(densum_local),COMM_S_NE)
      call mpireduce_sum_inplace(flowsum_local,shape(flowsum_local),COMM_S_NE) 
      ! misuse one line of the den_global array - its data is not
      ! needed any more
      call gather_array(den_global(:,1), number_of_species,&
         & densum_local, nsp, COMM_SP_NE)
      call gather_array(flow_global(:,1), number_of_species,&
         & flowsum_local, nsp, COMM_SP_NE)
    end if

    
    !Reduced over all directions 
    call mpiallreduce_sum_inplace(fluxsum1,nsp,COMM_SP_NE_S_NE)
    call mpiallreduce_sum_inplace(fluxsum2,1,COMM_X_EQ)
    call mpiallreduce_sum_inplace(fluxsum3,1,COMM_X_EQ)
    call mpiallreduce_sum_inplace(fluxsum2grad,1,COMM_X_EQ)
    !This one is only reduced along s as there is no velocity space
    !dependence.
    call mpiallreduce_sum_inplace(rmean,1,COMM_SP_NE_S_NE)
    call mpiallreduce_sum_inplace(phisum,5,COMM_S_NE)

    !Fourth order central difference
    gradphi = 8*(sqrt(phisum(4))-sqrt(phisum(2)))
    gradphi = gradphi + (-sqrt(phisum(5))+sqrt(phisum(1)))
    gradphi = gradphi*rrefoa/(12.E0*(neo_eps(2)-neo_eps(1)))

    !Only output the middle flux tube
    if(root_processor)write(file_unit,*) 'Flux average parflow <Bu_{||}>', flow_global(:,1) 
    if(root_processor)write(file_unit,*) 'Flux average density', den_global(:,1)
    if(root_processor)write(file_unit,*) 'Flux average potential', sqrt(phisum(3))
    if(root_processor)write(file_unit,*) 'Potential radial gradient', gradphi
    if(root_processor)write(file_unit,*) 'Diamagnetic toroidal flow', fluxsum2/rmean, &
       & fluxsum3/rmean
    if(root_processor)write(file_unit,*) 'Radial flow gradient', fluxsum2grad/rmean
    if(root_processor)write(file_unit,*) 'Aspect ratio, R/a', rrefoa

    if(root_processor)then
      write(*,*) 'Flux average parflow', flow_global(:,1) 
      write(*,*) 'Flux average density', den_global(:,1)
      write(*,*) 'Flux average potential', sqrt(phisum(3))
      write(*,*) 'Potential radial gradient', gradphi
      write(*,*) 'Diamagnetic toroidal flow', fluxsum2/rmean,fluxsum3/rmean
      write(*,*) 'Radial flow gradient', fluxsum2grad/rmean
      write(*,*) 'R-squared', rmean
    endif

    if(root_processor) then
      close(file_unit)
    end if
    if(allocated(den_global)) deallocate(den_global)
    if(allocated(flow_global)) deallocate(flow_global)

  end subroutine neof_diagnose

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_reconstr_distr_pot
    use io, only : output_array, xy_fmt, ascii_fmt
    use neoequil, only : neophi, neof
    use mpiinterface, only : root_processor

    if(root_processor) then
      ! FIXME gather these properly
      call output_array('neophi', 'diagnostic/diagnos_neoequil', &
         & neophi, 'F', &
         & xy_fmt, ascii_fmt)

      call output_array('neof', 'diagnostic/diagnos_neoequil', &
         & neof, 'F', &
         & xy_fmt, ascii_fmt)
    end if

  end subroutine output_reconstr_distr_pot

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_gradients
    use io, only : open_real_lu, close_lu, append_chunk
    use io, only : output_array, xy_fmt, ascii_fmt, attach_metadata
    use io, only : phys_unit_key, description_key, comments_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use neoequil, only : gradneophi, gradneof, neo_n_species_incl_dup
    use mpiinterface, only : root_processor, MPIREAL_X
    use mpiinterface, only : gather_array
    use mpidatatypes, only : create_subarray_datatype, free_datatype
    use global, only : id_s, id_mu, id_vpar
    use grid, only : proc_subset, ns, lsp, nsp
    use grid, only : n_s_grid, n_mu_grid, n_vpar_grid
    use mpicomms, only : COMM_S_NE, COMM_SP_EQ_X_EQ

    real, allocatable :: gradneof_global(:,:,:), gradneophi_global(:)
    integer :: mpi_dtype
    logical, parameter :: to_root_of_commcart = .true.
    logical :: root_is_in_slice

    integer :: lun_gradneof
    integer :: isg, dummy_isg

    if(root_processor) then

      call open_real_lu('gradneof', 'diagnostic/diagnos_neoequil', &
         & (/ n_s_grid, n_mu_grid, n_vpar_grid /), ascii_fmt, lun_gradneof)
      call attach_metadata_grid(lun_gradneof, 'neo_n_species_incl_dup', &
         & 's_grid', 'n_mu_grid', 'n_vpar_grid' , ascii_fmt)
      call attach_metadata(lun_gradneof, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(lun_gradneof, description_key, not_avail, ascii_fmt)
      call attach_metadata(lun_gradneof, comments_key, not_avail, ascii_fmt)

      
      allocate(gradneof_global(n_s_grid, n_mu_grid, &
         & n_vpar_grid))
      allocate(gradneophi_global(n_s_grid))
    else
      allocate(gradneof_global(1, 1, 1))
      allocate(gradneophi_global(1))
    end if
    if(proc_subset(1,0,0,0,0)) then
      call create_subarray_datatype(MPIREAL_X, mpi_dtype, &
         & id_s, id_mu, id_vpar)
      ! globally iterate over all species which have NEO backgrounds,
      ! and append one species after the other
      do isg = 1, neo_n_species_incl_dup
        if(proc_subset(1,0,0,0,isg) .or. root_processor) then
          root_is_in_slice = isg <= nsp
          ! root needs to pass a species index too, but lsp(isg) might
          ! yield sth < 1
          if(root_processor .and..not. root_is_in_slice) then
            dummy_isg = 1
          else
            dummy_isg = lsp(isg)
          end if
          call gather_array(gradneof_global, &
             & gradneof(dummy_isg,:,:,:), mpi_dtype, &
             & COMM_SP_EQ_X_EQ, to_root_of_commcart, root_is_in_slice)
          if(root_processor) then
            call append_chunk(lun_gradneof, gradneof_global, xy_fmt, ascii_fmt)
          end if
        end if
      end do
      call free_datatype(mpi_dtype)
    end if
    if(proc_subset(1,0,1,1,1)) then
      call gather_array(gradneophi_global, n_s_grid, &
         gradneophi, ns, COMM_S_NE)
    end if
    
    if(root_processor) then
      call output_array('gradneophi', 'diagnostic/diagnos_neoequil', &
         & gradneophi_global, 'F', &
         & xy_fmt, ascii_fmt)
    end if

    deallocate(gradneof_global)
    deallocate(gradneophi_global)

    if(root_processor) call close_lu(lun_gradneof, ascii_fmt)

  end subroutine output_gradients

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine output_neof_interp()
    use io, only : output_array, xy_fmt, ascii_fmt
    use neoequil, only : neo_coeffs_inter
    use mpiinterface, only : root_processor

    if(root_processor) then
      call output_array('neo_coeffs_inter', 'diagnostic/diagnos_neoequil', &
         & neo_coeffs_inter, 'F', &
         & xy_fmt, ascii_fmt)
    end if
    
  end subroutine output_neof_interp
  

  !--------------------------------------------------------------------
  !> If useful, print something to the stdout.
  !--------------------------------------------------------------------
  subroutine screen_output()
    use mpiinterface, only : root_processor

    if(.not. lcalc_neoequil) return
    
    if( root_processor) then
      !write(*,*)'Quantity X:', X
    end if
  end subroutine screen_output


end module diagnos_neoequil
