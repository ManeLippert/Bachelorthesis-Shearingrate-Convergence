!>------------------------------------------------------------------------------
!> ZONAL EVO MODULE DESCRIPTION
!> 
!>----DESCRIPTION---------------------------------------------------------------
!> 
!>   This diagnostic computes the time derivative of the zonal (flux surface
!>   averaged) electrostatic potential and the individual terms of the zonal
!>   evolution equation. 
!>
!>----OUTPUT--------------------------------------------------------------------
!>
!>   This diagnostic outputs several quantities connected to the evolution of
!>   the zonal (flux surface averaged) electrostatic potential.
!>
!>   The following output files have dimensions (time, kx) and
!>   contain the complex Fourier amplitudes of the respective quantity.
!> 
!>   If lcalc_zonal_evo = .true.
!>   - the zonal potential: zphi_kx
!>   - the time derivative of the zonal potential: dt_zphi_kx
!>   - all linear terms added: zevo_lin_kx
!>   - the nonlinear terms: zevo03_[phi,apar,bpar]_kx
!>
!>   If zevo_detail = .true.
!>   - the numerical damping terms: zevo_disp_[par,vp,perp]_kx
!>   - the contributions from linear terms: 
!>     zevo[01,02,04,05,07,08,10,11]_kx
!>   - the contribution from collisions: zevo_coll_kx
!>
!>   The following output files have dimensions (time, kx, ky) and
!>   contain the complex Fourier amplitudes of the respective quantity.
!> 
!>   If zevo_xy_spec = .true.
!>   - the nonlinear terms: zevo03_[phi,apar,bpar]_kxky
!>    
!>----FURTHER NOTES------------------------------------------------------------- 
!>
!>   Per default only the sum of all linear terms is calculated.
!>   More details about individual linear terms is obtaind by setting
!>   
!>       zevo_detail = .true.
!>
!>   kxky-spectra of nonlinear terms are output for
!> 
!>       zevo_xy_spec = .true.
!>
!-------------------------------------------------------------------------------
module diagnos_zonal_evo
  use control, only : lcalc_zonal_evo, zevo_detail
  
  implicit none

  private
  
  public :: set_default_nml_values, init, bcast, check, allocate_mem
  public :: finalize
  public :: initial_output
  public :: calc_smallstep
  public :: output
  public :: lcalc_zonal_evo, zevo_detail

  
  !> Outputs kxky-spectra of nonlinear terms
  logical, save, public :: zevo_xy_spec
  
  
  !> The logical unit number, used to output data
  integer, save :: i_zphi_kx = -1, i_dt_zphi_kx = -1
  integer, save :: i_zevo_lin_kx = -1
  integer, save :: i_zevo01_kx = -1, i_zevo02_kx = -1, i_zevo04_kx = -1
  integer, save :: i_zevo05_kx = -1, i_zevo07_kx = -1, i_zevo08_kx = -1
  integer, save :: i_zevo10_kx = -1, i_zevo11_kx = -1
  integer, save :: i_zevo_dpar_kx = -1, i_zevo_dperp_kx = -1, i_zevo_dvp_kx = -1
  integer, save :: i_zevo_coll_kx = -1
  integer, save :: i_zevo03_phi_kx = -1, i_zevo03_apar_kx = -1
  integer, save :: i_zevo03_bpar_kx = -1
  integer, save :: i_zevo03_phi_kxky = -1, i_zevo03_apar_kxky = -1
  integer, save :: i_zevo03_bpar_kxky = -1


  !> Generic integer used to select the appropriate nonlinear term
  integer, save :: NL_PHI  = 1
  integer, save :: NL_APAR = 2
  integer, save :: NL_BPAR = 3
  
  
  !> helper array for flux surface average
  complex, allocatable :: dum_za(:)
  
  !> helper array for the evaluation of the nonlinear term using
  !> hermitian symmery
  complex, allocatable :: dum_h(:,:)
  
  !> complex array 
  complex, allocatable :: zphi_kx(:), dt_zphi_kx(:)
  complex, allocatable :: zevo_dperp_kx(:), zevo_dvp_kx(:)
  complex, allocatable :: zevo01_kx(:), zevo02_kx(:)
  complex, allocatable :: zevo04_kx(:), zevo05_kx(:)
  complex, allocatable :: zevo07_kx(:), zevo08_kx(:)
  complex, allocatable :: zevo10_kx(:), zevo11_kx(:)
  complex, allocatable :: zevo_lin_kx(:)
  complex, allocatable :: zevo_dpar_kx(:), zevo_coll_kx(:)
  complex, allocatable :: zevo03_phi_kx(:), zevo03_apar_kx(:)
  complex, allocatable :: zevo03_bpar_kx(:)


  !> complex arrays living on kx,ky
  complex, allocatable :: zevo03_phi_kxky(:,:), zevo03_apar_kxky(:,:)
  complex, allocatable :: zevo03_bpar_kxky(:,:)
  
  !> arrays holding quantities of last time step to estimate time derivative
  complex, save, allocatable :: last_zphi_kx(:)
  
  !> time between call of calc_smallstep() and output()
  real, save :: delta_time_zonal_evo = 1.
  
  !> dummy array that holds the contribution of the individual terms of 
  !> GK equation to dt_g
  complex, allocatable :: dt_g(:)
  
  ! arrays for precalculated quantities
  complex, save, allocatable :: ci1kxrh(:)
  complex, save, allocatable :: ci1krho(:)
  complex, save, allocatable :: a_phi(:,:,:,:,:)
  complex, save, allocatable :: b_phi(:,:,:,:,:)
  complex, save, allocatable :: a_apar(:,:,:,:,:)
  complex, save, allocatable :: b_apar(:,:,:,:,:)
  complex, save, allocatable :: a_bpar(:,:,:,:,:)
  complex, save, allocatable :: b_bpar(:,:,:,:,:)
  
  
  complex, save, allocatable :: a(:,:)
  complex, save, allocatable :: b(:,:)
  complex, save, allocatable :: c(:,:)
  complex, save, allocatable :: d(:,:)

  complex, save, allocatable :: aa(:,:)
  complex, save, allocatable :: bb(:,:)
  complex, save, allocatable :: cc(:,:)
  complex, save, allocatable :: dd(:,:)

  real, save, allocatable :: ar(:,:)
  real, save, allocatable :: br(:,:)
  real, save, allocatable :: cr(:,:)
  real, save, allocatable :: dr(:,:)
  
  ! integer array to copy back the values of the arrays
  ! integer, allocatable, save :: lincopy(:)
  integer, save, allocatable :: lincopy3(:,:,:,:)

contains


  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()

    lcalc_zonal_evo = .false.
    zevo_detail = .false.
    zevo_xy_spec = .false.

  end subroutine set_default_nml_values


  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(lcalc_zonal_evo,1)
    call mpibcast(zevo_detail,1)
    call mpibcast(zevo_xy_spec,1)
    
  end subroutine bcast


  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()
    use general,     only : gkw_warn
    use control,     only : flux_tube, spectral_radius
    use source_time, only : source_time_ampl
    use global,      only : r_tiny
    use mode,        only : mode_box

    if(.not. lcalc_zonal_evo) return
    
    if(.not. flux_tube) then
      call gkw_warn('diagnos_zonal_evo requires flux_tube = .true..&
      & Diagnostic is switched off.')
      lcalc_zonal_evo = .false.
    end if
    
    if(.not. spectral_radius) then
      call gkw_warn('diagnos_zonal_evo requires spectral_radius = .true..&
      & Diagnostic is switched off.')
      lcalc_zonal_evo = .false.
    end if
    
    if(.not. mode_box) then
      call gkw_warn('diagnos_zonal_evo requires mode_box = .true..&
      & Diagnostic is switched off.')
      lcalc_zonal_evo = .false.
    end if
    
    if(zevo_detail) then
      call gkw_warn('The detailed output of diagnos_zonal_evo (zevo_detail)&
      & is memory consuming. ')
    end if
    
    if(zevo_xy_spec) then
      call gkw_warn('The kxky spectral output of diagnos_zonal_evo&
      & (zevo_xy_spec) is computationally expensive. If not mandatory,&
      & run performance might be improved by setting zevo_xy_spec = .false.')
    end if
    
    if(abs(source_time_ampl) > r_tiny) then
      call gkw_warn('diagnos_zonal_evo not implemented for finite&
      & energy source yet.&
      & Diagnostic is switched off.')
      lcalc_zonal_evo = .false.
    end if
    
  end subroutine check


  !--------------------------------------------------------------------
  !> Initialize the diagnostic. This is the place to open logical
  !> units.
  !--------------------------------------------------------------------
  subroutine init(requirements)
  
    use mpiinterface,    only : root_processor
    use io,              only : open_complex_lu, ascii_fmt, attach_metadata
    use io,              only : description_key, comments_key
    use io,              only : phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use diagnos_generic, only : LOCAL_DATA, S_GHOSTCELLS
    use diagnos_generic, only : MU_GHOSTCELLS, VPAR_GHOSTCELLS
    use global,          only : int2char_zeros
    use global,          only : PHI_FIELD, APAR_FIELD, DISTRIBUTION
    use global,          only : PHI_GA_FIELD, APAR_GA_FIELD, BPAR_GA_FIELD
    use grid,            only : nx, number_of_species, proc_subset
    use grid,            only : nmu, ns, nx, nmod, nsp, nvpar, n_s_grid
    use control,         only : non_linear, nlapar, nlphi, nlbpar
    use control,         only : lcollisions, disp_vp
    use functions,       only : besselj0_gkw, mod_besselj1_gkw
    use constants,       only : ci1
    use mode,            only : krho, kxrh, ixzero
    use components,      only : vthrat, tmp, signz, rhostar
    use velocitygrid,    only : mugr
    use fft,             only : four2D_real, four1D_real, FFT_FORWARD, FFT_INVERSE
    use non_linear_terms, only : mrad, jinv
    use dist,             only : ifdis
    use index_function,   only : indx
    
    logical, intent(inout) :: requirements(:,:)
    
    integer :: is, j, i, ix, imod, jv, kt, ipar
    integer :: idx2
    character (len=26) ::  filename
    real :: b0, b1_mod

    if(.not. lcalc_zonal_evo) return

    ! local data of distribution is required for all terms involving f/g
    requirements(DISTRIBUTION,LOCAL_DATA) = .true.
    
    ! s-ghostcells of distribution required for term I (and disp_par>0)
    requirements(DISTRIBUTION,S_GHOSTCELLS) = .true.
    
    ! vpar-ghostcells of distribution required for term IV (and disp_vp>0)
    requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
    
    ! local data of phi_ga required for term VII and VIII
    requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    
    ! s-ghostcells of phi_ga required for term VII
    requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
    
    ! local data of phi required for nonlinear term III
    if(nlphi) then
      requirements(PHI_FIELD,LOCAL_DATA) = .true.
    end if
    
    ! local data of apar_ga required for term III
    if(nlapar) then
      requirements(APAR_GA_FIELD,LOCAL_DATA) = .true.
    end if
    
    ! local data of bpar_ga required for term III
    if(nlbpar) then
      requirements(BPAR_GA_FIELD,LOCAL_DATA) = .true.
    end if
    
    
    
    if(lcollisions) then
        requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
        requirements(DISTRIBUTION,MU_GHOSTCELLS) = .true.
    end if
    if(disp_vp > 0) requirements(DISTRIBUTION,VPAR_GHOSTCELLS) = .true.
    
    
    ! precalculate quantities for optimization of nonlinear contribution
    if(non_linear) then
    
      ! precalculated quantities: complex krho and kxrh
      ci1krho(:) = ci1*krho(:)
      ci1kxrh(:) = ci1*kxrh(:)
  
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

      end do; end do; end do; end do; end do

      ! Call the FFTs once to force FFTW plan setup 
      a(:,:) = (0.0,0.0)
      ar(:,:) = 0.0
      call four2D_real(ar(:,:),a(:,:),FFT_INVERSE)
      call four2D_real(ar(:,:),a(:,:),FFT_FORWARD)
      
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

    end if
  
    
    
    if(root_processor) then
      
      ! zonal potential
      call open_complex_lu('zphi_kx', 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zphi_kx)
      call attach_metadata_grid(i_zphi_kx, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zphi_kx, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_zphi_kx, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_zphi_kx, comments_key, not_avail, ascii_fmt)
      
      ! time derivative of zonal potential
      call open_complex_lu('dt_zphi_kx', 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_dt_zphi_kx)
      call attach_metadata_grid(i_dt_zphi_kx, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_dt_zphi_kx, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_dt_zphi_kx, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_dt_zphi_kx, comments_key, not_avail, ascii_fmt)
      
      ! linear terms combined
      filename="zevo_lin_kx"
      call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
       & (/ nx /), ascii_fmt, i_zevo_lin_kx)
      call attach_metadata_grid(i_zevo_lin_kx, 'time', 'kxrh', ascii_fmt)
      call attach_metadata(i_zevo_lin_kx, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(i_zevo_lin_kx, description_key, not_avail, ascii_fmt)
      call attach_metadata(i_zevo_lin_kx, comments_key, not_avail, ascii_fmt)
      
      
      ! individual linear terms
      if(zevo_detail) then
      
        filename="zevo01_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo01_kx)
        call attach_metadata_grid(i_zevo01_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo01_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo01_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo01_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo02_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo02_kx)
        call attach_metadata_grid(i_zevo02_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo02_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo02_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo02_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo04_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo04_kx)
        call attach_metadata_grid(i_zevo04_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo04_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo04_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo04_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo05_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo05_kx)
        call attach_metadata_grid(i_zevo05_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo05_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo05_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo05_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo07_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo07_kx)
        call attach_metadata_grid(i_zevo07_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo07_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo07_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo07_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo08_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo08_kx)
        call attach_metadata_grid(i_zevo08_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo08_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo08_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo08_kx, comments_key, not_avail, ascii_fmt)
        
        if(nlbpar) then
          filename="zevo10_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo10_kx)
          call attach_metadata_grid(i_zevo10_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo10_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo10_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo10_kx, comments_key, not_avail, ascii_fmt)
          
          filename="zevo11_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo11_kx)
          call attach_metadata_grid(i_zevo11_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo11_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo11_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo11_kx, comments_key, not_avail, ascii_fmt)
        end if
        
        filename="zevo_dpar_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo_dpar_kx)
        call attach_metadata_grid(i_zevo_dpar_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo_dpar_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dpar_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dpar_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo_dperp_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo_dperp_kx)
        call attach_metadata_grid(i_zevo_dperp_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo_dperp_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dperp_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dperp_kx, comments_key, not_avail, ascii_fmt)
        
        filename="zevo_dvp_kx"
        call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
         & (/ nx /), ascii_fmt, i_zevo_dvp_kx)
        call attach_metadata_grid(i_zevo_dvp_kx, 'time', 'kxrh', ascii_fmt)
        call attach_metadata(i_zevo_dvp_kx, phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dvp_kx, description_key, not_avail, ascii_fmt)
        call attach_metadata(i_zevo_dvp_kx, comments_key, not_avail, ascii_fmt)
        
        if(lcollisions) then
          filename="zevo_coll_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo_coll_kx)
          call attach_metadata_grid(i_zevo_coll_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo_coll_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo_coll_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo_coll_kx, comments_key, not_avail, ascii_fmt)
        end if
      end if !zevo_detail
      
      
      ! nonlinear terms
      if(non_linear) then
      
        ! electrostatic exb-nonlinearity
        if(nlphi) then
        
          filename="zevo03_phi_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo03_phi_kx)
          call attach_metadata_grid(i_zevo03_phi_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo03_phi_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_phi_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_phi_kx, comments_key, not_avail, ascii_fmt)
          
          if(zevo_xy_spec) then
          
            filename="zevo03_phi_kxky"
            call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
             & (/ nx, nmod /), ascii_fmt, i_zevo03_phi_kxky)
            call attach_metadata_grid(i_zevo03_phi_kxky, 'time', 'kxrh', 'krho', ascii_fmt)
            call attach_metadata(i_zevo03_phi_kxky, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_phi_kxky, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_phi_kxky, comments_key, not_avail, ascii_fmt)
          
          end if
          
        end if
        
        
        if(nlapar) then
        
          ! electromagnetic magnetic flutter nonlinearity
          filename="zevo03_apar_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo03_apar_kx)
          call attach_metadata_grid(i_zevo03_apar_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo03_apar_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_apar_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_apar_kx, comments_key, not_avail, ascii_fmt)
          
          ! kxky-spectra
          if(zevo_xy_spec) then
          
            ! Maxwell stress
            filename="zevo03_apar_kxky"
            call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
             & (/ nx, nmod /), ascii_fmt, i_zevo03_apar_kxky)
            call attach_metadata_grid(i_zevo03_apar_kxky, 'time', 'kxrh', 'krho', ascii_fmt)
            call attach_metadata(i_zevo03_apar_kxky, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_apar_kxky, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_apar_kxky, comments_key, not_avail, ascii_fmt)
            
          end if
          
        end if
        
        
        if(nlbpar) then
        
          ! magnetic field compression
          filename="zevo03_bpar_kx"
          call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
           & (/ nx /), ascii_fmt, i_zevo03_bpar_kx)
          call attach_metadata_grid(i_zevo03_bpar_kx, 'time', 'kxrh', ascii_fmt)
          call attach_metadata(i_zevo03_bpar_kx, phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_bpar_kx, description_key, not_avail, ascii_fmt)
          call attach_metadata(i_zevo03_bpar_kx, comments_key, not_avail, ascii_fmt)
          
          ! kxky-spectra
          if(zevo_xy_spec) then
          
            ! magnetic field compression
            filename="zevo03_bpar_kxky"
            call open_complex_lu(filename, 'diagnostic/diagnos_zonal_evo', &
             & (/ nx, nmod /), ascii_fmt, i_zevo03_bpar_kxky)
            call attach_metadata_grid(i_zevo03_bpar_kxky, 'time', 'kxrh', 'krho', ascii_fmt)
            call attach_metadata(i_zevo03_bpar_kxky, phys_unit_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_bpar_kxky, description_key, not_avail, ascii_fmt)
            call attach_metadata(i_zevo03_bpar_kxky, comments_key, not_avail, ascii_fmt)
            
          end if
          
        end if
      
      end if
      
    end if

  end subroutine init


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()
  
    use grid,             only : nx, number_of_species
    use general,          only : gkw_abort
    use dist,             only : nsolc
    use control,          only : non_linear, nlapar, nlphi, nlbpar
    use control,          only : lcollisions
    use grid,             only : nmod, ns, nmu, nsp, nvpar, n_s_grid
    use non_linear_terms, only : mphiw3, mrad, mphi
  
  
    integer :: ierr
    

    if(.not. lcalc_zonal_evo) return
    
    
    ! allocate arrays for nonlinear contribution
    if(non_linear) then
    
      allocate(ci1krho(nmod),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate ci1krho in diagnos_zonal_evo')
      end if

      allocate(ci1kxrh(nx),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate ci1krho in diagnos_zonal_evo')
      end if
    
      allocate(a_phi(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate a_phi in diagnos_zonal_evo')
      end if
      
      allocate(b_phi(nmod,nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate a_phi in diagnos_zonal_evo')
      end if
      
      if(nlapar) then
      
        allocate(a_apar(nmod,nx,ns,nmu,nsp),stat=ierr)
        if (ierr /= 0) then
          call gkw_abort('Could not allocate a_apar in diagnos_zonal_evo')
        end if
        
        allocate(b_apar(nmod,nx,ns,nmu,nsp),stat=ierr)
        if (ierr /= 0) then
          call gkw_abort('Could not allocate a_apar in diagnos_zonal_evo')
        end if
      
      end if
      
      
      if(nlbpar) then
      
        allocate(a_bpar(nmod,nx,ns,nmu,nsp),stat=ierr)
        if (ierr /= 0) then
          call gkw_abort('Could not allocate a_bpar in diagnos_zonal_evo')
        end if
        
        allocate(b_bpar(nmod,nx,ns,nmu,nsp),stat=ierr)
        if (ierr /= 0) then
          call gkw_abort('Could not allocate a_bpar in diagnos_zonal_evo')
        end if
      
      end if


      allocate(a(mphiw3,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate a in diagnos_zonal_evo')
      end if
      
      allocate(b(mphiw3,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate b in diagnos_zonal_evo')
      end if
      
      allocate(c(mphiw3,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate c in diagnos_zonal_evo')
      end if
      
      allocate(d(mphiw3,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate d in diagnos_zonal_evo')
      end if
      
      allocate(aa(mrad,nmod), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate aa in diagnos_zonal_evo')
      end if
      
      allocate(bb(mrad,nmod), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate bb in diagnos_zonal_evo')
      end if
      
      allocate(cc(mrad,nmod), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate cc in diagnos_zonal_evo')
      end if
      
      allocate(dd(mrad,nmod), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate dd in diagnos_zonal_evo')
      end if
      
      allocate(ar(mphi,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate ar in diagnos_zonal_evo')
      end if
      
      allocate(br(mphi,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate br in diagnos_zonal_evo')
      end if

      allocate(cr(mphi,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate cr in diagnos_zonal_evo')
      end if
      
      allocate(dr(mphi,mrad), stat = ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate dr in diagnos_zonal_evo')
      end if
      
      allocate(lincopy3(nvpar*nmod*nx,ns,nmu,nsp),stat=ierr)
      if (ierr /= 0) then
        call gkw_abort('Could not allocate lincopy in non_linear_terms')
      end if
          
    end if
    
    
    allocate(zphi_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zphi_kx')
    
    allocate(last_zphi_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array last_zphi_kx')
    
    allocate(dt_zphi_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array dt_zphi_kx')
    
    allocate(zevo_lin_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo_lin_kx')
  
    allocate(zevo01_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo01_kx')
  
    allocate(zevo02_kx(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo02_kx')
  
    if(zevo_detail) then
      allocate(zevo04_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo04_kx')
    
      allocate(zevo05_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo05_kx')
      
      allocate(zevo07_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo07_kx')
    
      allocate(zevo08_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo08_kx')
      
      if(nlbpar) then
        allocate(zevo10_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo10_kx')
        
        allocate(zevo11_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo11_kx')
      end if
      
      allocate(zevo_dpar_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo_dpar_kx')
      
      allocate(zevo_dvp_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo_dvp_kx')
      
      allocate(zevo_dperp_kx(nx),stat=ierr)
      if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo_dperp_kx')
      
      if(lcollisions) then
        allocate(zevo_coll_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo_coll_kx')
      end if
    end if ! zevo_detail
    
    
    allocate(dum_za(nx),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array dum_za')
    
    allocate(dum_h(mrad,nmod),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array dum_h')
    
    
    if(non_linear) then
    
      if(nlphi) then
        allocate(zevo03_phi_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_phi_kx')
        
        if(zevo_xy_spec) then
          allocate(zevo03_phi_kxky(nx,nmod),stat=ierr)
          if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_phi_kxky')
        end if
        
      end if
      
      if(nlapar) then
        allocate(zevo03_apar_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_apar_kx')
        
        if(zevo_xy_spec) then
          allocate(zevo03_apar_kxky(nx,nmod),stat=ierr)
          if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_apar_kxky')
        end if
        
      end if
      
      if(nlbpar) then
        allocate(zevo03_bpar_kx(nx),stat=ierr)
        if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_bpar_kx')
        
        if(zevo_xy_spec) then
          allocate(zevo03_bpar_kxky(nx,nmod),stat=ierr)
          if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array zevo03_bpar_kxky')
        end if
        
      end if
      
    end if
    
    allocate(dt_g(nsolc),stat=ierr)
    if (ierr /= 0) call gkw_abort('zonal_evo :: could not allocate array dt_g')

    
  end subroutine allocate_mem



  !--------------------------------------------------------------------
  !> Clean up, deallocate, close everything.
  !--------------------------------------------------------------------
  subroutine finalize()
  
    use io,            only : close_lu, ascii_fmt
    use grid,          only : number_of_species
    use mpiinterface,  only : root_processor
    use control,       only : non_linear, nlapar, nlphi, nlbpar
    use control,       only : lcollisions
    
    if(.not. lcalc_zonal_evo) return
    
    ! deallocate arrays if allocated
    if(allocated(zphi_kx)) deallocate(zphi_kx)
    if(allocated(dt_zphi_kx)) deallocate(dt_zphi_kx)
    
    if(allocated(zevo_lin_kx)) deallocate(zevo_lin_kx)
    if(allocated(zevo01_kx)) deallocate(zevo01_kx)
    if(allocated(zevo02_kx)) deallocate(zevo02_kx)
    if(allocated(zevo04_kx)) deallocate(zevo04_kx)
    if(allocated(zevo05_kx)) deallocate(zevo05_kx)
    if(allocated(zevo07_kx)) deallocate(zevo07_kx)
    if(allocated(zevo08_kx)) deallocate(zevo08_kx)
    if(allocated(zevo10_kx)) deallocate(zevo10_kx)
    if(allocated(zevo11_kx)) deallocate(zevo11_kx)
    if(allocated(zevo_dpar_kx)) deallocate(zevo_dpar_kx)
    if(allocated(zevo_dvp_kx)) deallocate(zevo_dvp_kx)
    if(allocated(zevo_dperp_kx)) deallocate(zevo_dperp_kx)
    if(allocated(zevo_coll_kx)) deallocate(zevo_coll_kx)
    if(allocated(zevo03_phi_kx)) deallocate(zevo03_phi_kx)
    if(allocated(zevo03_apar_kx)) deallocate(zevo03_apar_kx)
    if(allocated(zevo03_bpar_kx)) deallocate(zevo03_bpar_kx)
    
    if(allocated(zevo03_phi_kxky)) deallocate(zevo03_phi_kxky)
    if(allocated(zevo03_apar_kxky)) deallocate(zevo03_apar_kxky)
    if(allocated(zevo03_bpar_kxky)) deallocate(zevo03_bpar_kxky)
    
    if(allocated(dt_g)) deallocate(dt_g)
    
    if(allocated(dum_za)) deallocate(dum_za)
    if(allocated(dum_h)) deallocate(dum_h)
    
    if(allocated(a)) deallocate(a)
    if(allocated(b)) deallocate(b)
    if(allocated(c)) deallocate(c)
    if(allocated(d)) deallocate(d)
    if(allocated(aa)) deallocate(aa)
    if(allocated(bb)) deallocate(bb)
    if(allocated(cc)) deallocate(cc)
    if(allocated(dd)) deallocate(dd)
    if(allocated(ar)) deallocate(ar)
    if(allocated(br)) deallocate(br)
    if(allocated(cr)) deallocate(cr)
    if(allocated(dr)) deallocate(dr)
    
    if(allocated(a_phi)) deallocate(a_phi)
    if(allocated(b_phi)) deallocate(b_phi)
    if(allocated(a_apar)) deallocate(a_apar)
    if(allocated(b_apar)) deallocate(b_apar)
    if(allocated(a_bpar)) deallocate(a_bpar)
    if(allocated(b_bpar)) deallocate(b_bpar)
    
    if(allocated(ci1krho)) deallocate(ci1krho)
    if(allocated(ci1kxrh)) deallocate(ci1kxrh)
    
    if(allocated(lincopy3)) deallocate(lincopy3)
    
    
    ! Close all logical units
    if(root_processor) then
    
      call close_lu(i_zphi_kx, ascii_fmt)
      call close_lu(i_dt_zphi_kx, ascii_fmt)
      
      call close_lu(i_zevo_lin_kx, ascii_fmt)
        
      if(zevo_detail) then
        call close_lu(i_zevo01_kx, ascii_fmt)
        call close_lu(i_zevo02_kx, ascii_fmt)
        call close_lu(i_zevo04_kx, ascii_fmt)
        call close_lu(i_zevo05_kx, ascii_fmt)
        call close_lu(i_zevo07_kx, ascii_fmt)
        call close_lu(i_zevo08_kx, ascii_fmt)
        if(nlbpar) then
          call close_lu(i_zevo10_kx, ascii_fmt)
          call close_lu(i_zevo11_kx, ascii_fmt)
        end if
        call close_lu(i_zevo_dpar_kx, ascii_fmt)
        call close_lu(i_zevo_dperp_kx, ascii_fmt)
        call close_lu(i_zevo_dvp_kx, ascii_fmt)
        if(lcollisions) then
          call close_lu(i_zevo_coll_kx, ascii_fmt)
        end if
      end if
      
      if(non_linear) then
      
        if(nlphi) then
        
          call close_lu(i_zevo03_phi_kx, ascii_fmt)
          
          if (zevo_xy_spec) then
            call close_lu(i_zevo03_phi_kxky, ascii_fmt)
          end if
          
        end if
        
        if(nlapar) then
        
          call close_lu(i_zevo03_apar_kx, ascii_fmt)
          
          if(zevo_xy_spec) then
            call close_lu(i_zevo03_apar_kxky, ascii_fmt)
          end if
          
        end if
        
        if(nlbpar) then
        
          call close_lu(i_zevo03_bpar_kx, ascii_fmt)
          
          if(zevo_xy_spec) then
            call close_lu(i_zevo03_bpar_kxky, ascii_fmt)
          end if
          
        end if
      
      end if
    
    end if

    
  end subroutine finalize

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()

    if(.not. lcalc_zonal_evo) return

  end subroutine initial_output


  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine calc_smallstep(i_smallstep)
  
    use control,          only : naverage, time
    real, save :: last_time_pre_zonal_evo = 0.
    
    integer, intent(in) :: i_smallstep

    if (.not.lcalc_zonal_evo) return

    if (i_smallstep == naverage - 1) then
    
      ! calculate zonal phi at time step just before big time step, in order
      ! to estimate time derivative
      call calc_zphi(last_zphi_kx)
      ! store the time this subroutine was last called in a local variable
      last_time_pre_zonal_evo = time

    else if (i_smallstep == naverage) then
      ! save time interval between this step and the
      ! last call of calc_smallstep().
      delta_time_zonal_evo = time - last_time_pre_zonal_evo
    end if
    
  end subroutine calc_smallstep


  !--------------------------------------------------------------------
  !> The routine calc_largestep() is to be called repeatedly, after every
  !> large timestep, and here the diagnostic should calculate its
  !> quantities.
  !> 
  !> Splitting calculation and output is to be encouraged, particularly when
  !>  - the output is higher dimensional, and requires MPI or MPI-IO
  !>  - same data is output in different dimensional slices
  !--------------------------------------------------------------------
  subroutine calc_largestep()
    use grid,           only : proc_subset
    use grid,           only : number_of_species, nx, lsp
    use mpiinterface,   only : mpireduce_sum_inplace
    use global,         only : r_tiny
    use constants,      only : c1
    use dist,           only : ifdis, iphi
    use matdat,         only : matd, matvpd, matperpd, matcoll
    use matdat,         only : mat_vpar_grd_phi, mat_vd_grad_phi_fm
    use matdat,         only : mat_vpar_grad_df, mat_vdgradf
    use matdat,         only : mat_trapdf_4d, mat_ve_grad_fm
    use matdat,         only : mat_vpar_grd_phi, mat_vd_grad_phi_fm
    use matdat,         only : mat_vpar_grd_bpar, mat_vd_grad_bpar_fm
    use matdat,         only : mat, mat_maxwll_background
    use index_function, only : indx
    use general,        only : gkw_abort
    use matrix_format,  only : usmv
    use components,     only : adiabatic_electrons
    use control,        only : zonal_adiabatic, non_linear, nlapar
    use control,        only : nlphi, nlbpar, lcollisions
    use fields,         only : calculate_fields
    
    
    integer :: ix
    
    if(.not. lcalc_zonal_evo) return
    
    ! calculate zonal potential spectrum
    call calc_zphi(zphi_kx)
    
    ! calculate the time derivative of the zonal potential
    if(proc_subset(0,1,1,1,1)) then
      do ix=1, nx
        if (abs(zphi_kx(ix) - last_zphi_kx(ix)) < r_tiny) then
          dt_zphi_kx(ix) = 0.0
        else
          dt_zphi_kx(ix) = (zphi_kx(ix) - &
               & last_zphi_kx(ix))/delta_time_zonal_evo
        end if
      end do
    end if

    ! Calculate individual linear contributions.
    ! Notation from the GKW documentation is applied: 
    !   'I : vpar_grad_df'
    !   'II : vdgradf'
    !   'IV : trapdf_4d'
    !   'V : ve_grad_fm'
    !   'VII : vpar_grd_phi (Landau damping)'
    !   'VIII : vd_grad_phi_fm'
    !   'X : Trapping due to the perturbed magnetic field'
    !   'XI : vd_grad_bpar_fm'
    
    if(zevo_detail) then
      ! term I
      call calc_linear_contribution(zevo01_kx, mat_vpar_grad_df)
      
      ! term II
      call calc_linear_contribution(zevo02_kx, mat_vdgradf)
      
      ! term IV
      call calc_linear_contribution(zevo04_kx, mat_trapdf_4d)
      
      ! term V
      call calc_linear_contribution(zevo05_kx, mat_ve_grad_fm)
      
      ! term VII
      call calc_linear_contribution(zevo07_kx, mat_vpar_grd_phi)
      
      ! term VIII
      call calc_linear_contribution(zevo08_kx, mat_vd_grad_phi_fm)
      
      if(nlbpar) then
        ! term X
        call calc_linear_contribution(zevo10_kx, mat_vpar_grd_bpar)
        
        ! term XI
        call calc_linear_contribution(zevo11_kx, mat_vd_grad_bpar_fm)
      end if
      
      if(lcollisions) then
        ! collisions
        call calc_linear_contribution(zevo_coll_kx, matcoll)
      end if
      
      ! disp_par
      call calc_linear_contribution(zevo_dpar_kx, matd)
      
      ! disp_perp
      call calc_linear_contribution(zevo_dperp_kx, matperpd)
      
      ! disp_vp
      call calc_linear_contribution(zevo_dvp_kx, matvpd)
    
      ! add all linear terms
      zevo_lin_kx = (0.E0, 0.E0)
      do ix=1, nx
        zevo_lin_kx(ix) = zevo01_kx(ix)+zevo02_kx(ix)+zevo04_kx(ix) &
        & + zevo05_kx(ix)+zevo07_kx(ix)+zevo08_kx(ix) &
        & + zevo_dpar_kx(ix)+zevo_dperp_kx(ix)+zevo_dvp_kx(ix)
        
        if(lcollisions) then
          zevo_lin_kx(ix) = zevo_lin_kx(ix) + zevo_coll_kx(ix)
        end if
        
        if(nlbpar) then
          zevo_lin_kx(ix) = zevo_lin_kx(ix) + zevo10_kx(ix) + zevo11_kx(ix)
        end if
      end do
      
    else ! zevo_detail
    
      ! Use zevo01_kx as dummy for all linear terms involving the 
      ! distribution function.
      call calc_linear_contribution(zevo01_kx, mat)
      
      ! Use zevo02_kx as dummy for all linear terms involving the fields.
      call calc_linear_contribution(zevo02_kx, mat_maxwll_background)
      
      ! add all linear terms
      zevo_lin_kx = (0.E0, 0.E0)
      do ix=1, nx
        zevo_lin_kx(ix) = zevo01_kx(ix)+zevo02_kx(ix)
      end do
    
    end if ! .not. zevo_detail
    
    
    ! Caluclate nonlinear contributions
    if(non_linear) then
    
      if(nlphi) then
        ! electrostatic exb-nonlinearity
        call calc_nonlinear_contribution(zevo03_phi_kx, NL_PHI)
        
        if(zevo_xy_spec) then
          call calc_nonlin_ctrb_hermitian(zevo03_phi_kxky, NL_PHI)
        end if
        
      end if
      
      if(nlapar) then
        ! electromagnetic magnetic flutter nonlinearity
        call calc_nonlinear_contribution(zevo03_apar_kx, NL_APAR)
        
        if(zevo_xy_spec) then
          call calc_nonlin_ctrb_hermitian(zevo03_apar_kxky, NL_APAR)
        end if
        
      end if
      
      if(nlbpar) then
        ! electromagnetic correction due to parallel magnetic field
        ! compression
        call calc_nonlinear_contribution(zevo03_bpar_kx, NL_BPAR)
        
        if(zevo_xy_spec) then
          call calc_nonlin_ctrb_hermitian(zevo03_bpar_kxky, NL_BPAR)
        end if
        
      end if
    
    end if
    
  end subroutine calc_largestep


  !--------------------------------------------------------------------
  !> Subroutine that calculates the radial spectrum of the zonal
  !> electrostatic potential
  !--------------------------------------------------------------------
  subroutine calc_zphi(zphi)
  
    use dist,           only : fdisi, iphi
    use grid,           only : nx, ns, proc_subset
    use geom,           only : ints
    use mode,           only : iyzero
    use mpiinterface,   only : mpireduce_sum_inplace
    use mpicomms,       only : COMM_S_NE
    use index_function, only : indx
    
    complex, intent(out) :: zphi(:)
    integer :: i, ix
    
    
    ! calculate zonal (flux surface averaged) potential
    if(proc_subset(0,0,1,1,1)) then
      zphi = 0.0
      do ix = 1, nx
        do i = 1, ns
            zphi(ix) = zphi(ix) + ints(i) * fdisi(indx(iphi,iyzero,ix,i))
        end do
      end do
      ! complete flux surface average
      call mpireduce_sum_inplace(zphi,shape(zphi),COMM_S_NE)
    end if
    
  end subroutine calc_zphi



  
  !--------------------------------------------------------------------
  !> This subroutine calculates the contribution of a linear term to 
  !> the evolution of the zonal potential as follows:
  !> 
  !> (i) Calculate contribution to dg/dt via unstructured sparse matrix
  !>     vector multiplication (usmv)
  !> (ii) Inverts Poisson equation with dg/dt replacing distribution
  !>      function in integral part  (calculate_fields)
  !> (iii) Calculate the flux surface average
  !--------------------------------------------------------------------
  subroutine calc_linear_contribution(ctrb, matr)
  
    use grid,           only : nx, ns, proc_subset
    use geom,           only : ints
    use mode,           only : iyzero
    use fields,         only : calculate_fields
    use dist,           only : iphi, fdis_tmp
    use index_function, only : indx
    use mpicomms,       only : COMM_S_NE
    use mpiinterface,   only : mpireduce_sum_inplace
    use matrix_format,  only : usmv, sparse_matrix
    use constants,      only : c1
    
    
    complex, intent(out) :: ctrb(:)
    type(sparse_matrix), intent(in) :: matr
    
    integer :: ix, i, iih, ierr
    
        
    ! calculate dt_g for given linear term (fdis_tmp contains f not g
    ! and, hence, its usage is correct for linear terms.)
    dt_g = 0.0
    call usmv(c1,matr,fdis_tmp,dt_g,ierr)
    
  
    ! invert poisson equation using dt_g as distribution in its 
    ! integral part
    call calculate_fields(dt_g)
    
    ! Calculate flux-surface average
    ! Work on spatial processor subset only, since electrostatic potential
    ! has no velocity space- and species dependence.
    if(proc_subset(0,0,1,1,1)) then
    
      ! perform flux surface average
      ctrb = 0.0
      do ix=1, nx
        do i=1, ns
        
          ! index of phi
          iih = indx(iphi,iyzero,ix,i) 
          
          ! integral in s
          ctrb(ix) = ctrb(ix) + ints(i) * dt_g(iih)
          
        end do
      end do
      
      !finish flux surface average
      call mpireduce_sum_inplace(ctrb,shape(ctrb), COMM_S_NE)
    
    end if
  
  end subroutine calc_linear_contribution



  !-----------------------------------------------------------------------------
  !> This Subroutine computes the nonlinear contribution, i.e., the term 
  !> v_\chi \dot \nabla with v_\chi, to the zonal potential evolution.
  !> Calculation is splitted into electrostatic and electromagnetic parts.
  !> The respective part is selected by
  !>            stress = NL_PHI -> v_\chi -> exb-drift
  !>            stress = NL_APAR -> v_\chi -> parallel motion along
  !>                                        pertrubed magnetic field lines,
  !>                                        i.e., magnetic flutter 
  !>            stress = NL_BPAR -> v_\chi -> magnetic field compression
  !> Calculation is done as follows:
  !> 
  !> (i) Calculate nonlinear contribution to dg/dt via semi-spectral 
  !>     methods.
  !> (ii) Invert Poisson equation with dg/dt replacing distribution
  !>      function in integral part  (calculate_fields)
  !> (iii) Calculate the flux surface average.
  !-----------------------------------------------------------------------------
  subroutine calc_nonlinear_contribution(ctrb,stress)

    use dist,             only : get_phi, phi, get_apar, apar
    use dist,             only : fdisi, ifdis, iphi, iapar, ibpar
    use dist,             only : get_bpar, bpar
    use mode,             only : ixzero, iyzero
    use geom,             only : efun, ints
    use control,          only : non_linear, nlapar, nlphi, nlbpar
    use grid,             only : nmod, nx, ns, nmu, nvpar, nsp, proc_subset
    use velocitygrid,     only : vpgr
    use fft,              only : four2D_real, FFT_INVERSE, FFT_FORWARD
    use general,          only : gkw_abort
    use global,           only : r_tiny, r_huge
    use constants,        only : c1
    use non_linear_terms, only : jind, mphi, mrad
    use index_function,   only : indx
    use mpicomms,         only : COMM_S_NE
    use mpiinterface,     only : mpireduce_sum_inplace
    use fields,           only : calculate_fields

    complex, intent(out) :: ctrb(:)
    integer, intent(in) :: stress

    complex :: cdum
    integer :: idx, idxcopy, idx2, idxcopy2
    integer :: imod, ix, i, j, is, iih
    integer :: jv, ipar, kt
    real    :: dum, mphimrad1

    
    ! (i) Calculate nonlinear contribution to dg/dt via semi-spectral
    !     methods. 
    
    ! use mphimrad1i =  c1 / (mrad*mphi)
    mphimrad1 = c1 / (mrad*mphi)
    

    ! Abort if called incorrectly
    if (.not. (non_linear)) then
      call gkw_abort('Invalid call to calc_nonlinear_contribution')
    end if


    ! update the fields in fdisi according to nonlinear contribution
    
    ! electrostatic or Reynolds stress
    if (nlphi .and. stress == NL_PHI) then
      call get_phi(fdisi,phi)
    end if
    
    ! electromagnetic or Maxwell stress
    if (nlapar .and. stress == NL_APAR) then
      call get_apar(fdisi,apar)
    end if
    
    ! correction due to parallel field compression
    if (nlbpar .and. stress == NL_BPAR) then
      call get_bpar(fdisi,bpar)
    end if 

    ! initialize the indices of the help arrays
    idx = 1
    idxcopy = 1
    
    ! set dummy for dt_g to zero
    dt_g = (0.0, 0.0)
    
    do ipar=1, ns; do is=1, nsp; do jv=1, nmu
    
      ! Gyroaveraged field in k space obtained with Bessel function J_0
      ! <field_k> = J_0(k_perp rho) field_k
      ! Does not depend on parallel velocity, calculated outside loop

      ! a=grad_y_k <field_k> = zeta gradient of the gyroaverage of field
      !                      = i J_0() k_zeta field_k
      a(:,:) = (0.,0.)

      ! b=grad_x_k <field_k> = radial gradient of the gyroaverage of field
      !                     = i J_0() k_psi field_k
      b(:,:) = (0.,0.)
      
      ! below field is chosen among electrostatic potential phi (NL_PHI)
      ! or parallel vector potential apar (NL_APAR)
      
      if(stress == NL_PHI) then

        loop_nx1: do ix = 1, nx
          loop_nmod1: do imod = 1, nmod
            a(imod,jind(ix)) = a_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
            b(imod,jind(ix)) = b_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
          end do loop_nmod1
        end do  loop_nx1
      
      end if
      
      
      if(stress == NL_APAR) then

        do ix = 1, nx        
          do imod = 1, nmod  
            a(imod,jind(ix)) = a_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
            b(imod,jind(ix)) = b_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
          end do
        end do

      end if
      
      
      if(stress == NL_BPAR) then

        do ix = 1, nx        
          do imod = 1, nmod  
            a(imod,jind(ix)) = a_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
            b(imod,jind(ix)) = b_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
          end do
        end do

      end if


      ! Inverse fourier transform of field from k space to real space
      ! ar = grad_zeta <field> =
      ! poloidal gradient of the gyroaverage of field in real space
      ! br = grad_psi <field> =
      ! radial gradient of the gyroaverage of field in real space
      call four2D_real(ar(:,:),a(:,:),FFT_INVERSE)
      call four2D_real(br(:,:),b(:,:),FFT_INVERSE)
    

      ! Not until here is the parallel velocity looped over
      idxcopy2 = 1
      idx2 = 1

      loop_nvpar: do kt = 1, nvpar

        ! initialize a and b
        a(:,:) = (0.,0.)
        b(:,:) = (0.,0.)

        ! a = grad_x f_k
        ! b = grad_y f_k
        loop_nx2: do ix = 1, nx            ! Loop over radial modes
          cdum = ci1kxrh(ix)
          loop_nmod2: do imod = 1, nmod    ! Loop over poloidal modes

            a(imod,jind(ix)) =  cdum          * fdisi(indx(ifdis,imod,ix,ipar,jv,kt,is))
            b(imod,jind(ix)) =  ci1krho(imod) * fdisi(indx(ifdis,imod,ix,ipar,jv,kt,is))
            idx2 = idx2 + 1
          end do loop_nmod2
        end do loop_nx2

        ! Inverse fourier transform from k space to real space
        ! For gradients of distribution function
        call four2D_real(cr(:,:),a(:,:),FFT_INVERSE)
        call four2D_real(dr(:,:),b(:,:),FFT_INVERSE)

        ! Now everything is in real space
        ! cr = grad_x f = radial gradient of the distribution in real space
        ! dr = grad_y f = zeta gradient of the distribution in real space

        ! The zeta psi component = efun(ipar,2,1)
        ! The psi zeta component = -efun(ipar,2,1) = efun(ipar,1,2)
        ! Use the efun of the first radial grid point, since no radial
        ! dependence in flux-tube version.
        if(stress == NL_PHI) then
          dum = -efun(1,ipar,2,1)
        end if
        ! No minus sign, since apar correction is negative
        if(stress == NL_APAR) then
          dum = vpgr(ipar,jv,kt,is) * efun(1,ipar,2,1)
        end if
        if(stress == NL_BPAR) then
          dum = -efun(1,ipar,2,1)
        end if

        ! Evaluate nonlinear term in real space
        do j = 1, mrad
          do i = 1, mphi
            ! The minus sign is due to the antisymmetry of efun
            ! efun is antisymmetric for all geometries, not just circular
            cr(i,j) = dum*(ar(i,j)*cr(i,j)-br(i,j)*dr(i,j))
          end do
        end do


        ! Forward Fourier transform from real back to k space.
        ! Note that the positions of input and output arguments are reversed.
        call four2D_real(cr(:,:),a(:,:),FFT_FORWARD)
        
        ! Normalise
        a(:,:) = mphimrad1*a(:,:)
        ! Now the array a contains the chosen nonlinear term,
        ! evaluated for this timestep and transformed to k space.
  
      
        ! Copy the result of the calculation into the dt_g array:
        ! First start the i index from 2 to avoid adding in the (0,0) mode
        j = 1
        idxcopy2 = idxcopy2 + 1 ! start at 1 to avoid the first element
        do i = 2, nmod
          dt_g(lincopy3(idxcopy2,ipar,jv,is)) = dt_g(lincopy3(idxcopy2,ipar,jv,is)) + a(i,j)
          idxcopy2 = idxcopy2 + 1
        end do

        do j = 2, nx-ixzero+1
          do i = 1, nmod
            dt_g(lincopy3(idxcopy2,ipar,jv,is)) = dt_g(lincopy3(idxcopy2,ipar,jv,is)) + a(i,j)
            idxcopy2 = idxcopy2 + 1
          end do
        end do

        do j = mrad+2-ixzero, mrad
          do i = 1, nmod
            dt_g(lincopy3(idxcopy2,ipar,jv,is)) = dt_g(lincopy3(idxcopy2,ipar,jv,is)) + a(i,j)
            idxcopy2 = idxcopy2 + 1
          end do
        end do

      end do loop_nvpar

    end do; end do; end do 
    
    
    !(ii) Invert Poisson equation with dg/dt replacing distribution
    !      function in integral part.
    call calculate_fields(dt_g)
    
   
    ! (iii) Calculate the flux surface average.
    if(proc_subset(0,0,1,1,1)) then
      ctrb = (0.E0, 0.E0)
      do ix=1, nx
        do i=1, ns
        
          ! index of phi
          iih = indx(iphi,iyzero,ix,i) 
          
          ! integral in s
          ctrb(ix) = ctrb(ix) + ints(i) * dt_g(iih)
          
        end do
      end do
      
      !finish flux surface average
      call mpireduce_sum_inplace(ctrb,shape(ctrb), COMM_S_NE)
    end if
   
  end subroutine calc_nonlinear_contribution
  
  

  
  
  !-----------------------------------------------------------------------------
  !> This Subroutine computes the nonlinear contribution, i.e. the term 
  !> v_\chi \dot \nabla with v_\chi, to the zonal potential evolution.
  !> Calculation is splitted into electrostatic and electromagnetic part.
  !> The respective part is selected by
  !>            stress = NL_PHI -> v_\chi -> exb-drift
  !>            stress = NL_APAR -> v_\chi -> parallel motion along
  !>                                        pertrubed magnetic field lines,
  !>                                        i.e., magnetic flutter 
  !>            stress = NL_BPAR -> v_\chi -> magnetic field compression
  !> This routine uses the Hermitian symmetry of the spectrum together with
  !> the flux-surface average in order to decompose the nonlinear contribution
  !> at fixed radial mode into binormal modes.
  !> Calculation is done as follows:
  !> 
  !> (i) Calculate nonlinear contribution to dg/dt via semi-spectral 
  !>     methods. Fourier transform in radial direction only, since 
  !>     number of connecting kzeta-modes collapses due to flux-surface
  !>     average over squared quantity.
  !> (ii) Invert Poisson equation with dg/dt replacing distribution
  !>      function in integral part  (calculate_fields)
  !> (iii) Calculate the flux surface average.
  !-----------------------------------------------------------------------------
  subroutine calc_nonlin_ctrb_hermitian(ctrb,stress)

    use dist,             only : get_phi, phi, get_apar, apar
    use dist,             only : ifdis, iphi, iapar, fdisi, ibpar
    use dist,             only : get_bpar, bpar
    use mode,             only : ixzero, iyzero
    use geom,             only : efun, ints
    use control,          only : non_linear, nlapar, nlphi
    use control,          only : nlbpar, spectral_radius
    use grid,             only : nmod, nx, ns, nmu, nvpar, nsp, proc_subset
    use velocitygrid,     only : vpgr
    use fft,              only : FFT_INVERSE, FFT_FORWARD, fourcol
    use general,          only : gkw_abort
    use constants,        only : c1
    use mpiinterface,     only : mpireduce_sum_inplace
    use functions,        only : besselj0_gkw
    use non_linear_terms, only : jind, jinv, mrad
    use index_function,   only : indx
    use mpicomms,         only : COMM_S_NE
    use fields,           only : calculate_fields
    use diagnos_generic,  only : parseval_correction


    ! Here ctrb has two dimensions (nx,nmod)
    complex, intent(out) :: ctrb(:,:)
    integer, intent(in) :: stress

    integer :: imod, ix, i, j, is, iih
    integer :: jv, ipar, kt
    real    :: dum


    ! Abort if called incorrectly
    if (.not. (non_linear)) then
      call gkw_abort('Invalid call to calc_nonlinear_contribution')
    end if
    
    
    ! (i) Calculate nonlinear contribution to dg/dt via semi-spectral 
    !     methods.
    
    ! update the fields in the phi, apar according to stress
    ! electrostatic or Reynolds stress
    if (nlphi .and. stress == NL_PHI) then
      call get_phi(fdisi,phi)
    end if
    
    ! electromagnetic or Maxwell stress
    if (nlapar .and. stress == NL_APAR) then
      call get_apar(fdisi,apar)
    end if
    
    ! correction due to parallel field compression
    if (nlbpar .and. stress == NL_BPAR) then
      call get_bpar(fdisi,bpar)
    end if
    
    
    ! Calculate contribution per kzeta mode
    loop_nmod: do imod=1, nmod
     
      ! set dummy for dt_g to zero
      dt_g = (0.0, 0.0)
      
      do ipar=1, ns; do is=1, nsp; do jv=1, nmu
      
        ! Gyroaveraged field in k space obtained with Bessel function J_0
        ! <field_k> = J_0(k_perp rho) field_k
        ! Does not depend on parallel velocity, calculated outside loop

        ! aa=grad_y_k <field_k> = zeta gradient of the gyroaverage of field
        !                      = i J_0() k_zeta field_k
        aa(:,:) = (0.,0.)

        ! bb=grad_x_k <field_k> = radial gradient of the gyroaverage of field
        !                     = i J_0() k_psi field_k
        bb(:,:) = (0.,0.)
        
        
        ! below field is chosen among electrostatic potential phi (NL_PHI)
        ! or parallel vector potential apar (NL_APAR)
        loop_nx1: do ix=1, nx
        
          if(stress == NL_PHI) then
            aa(jind(ix),imod) = a_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
            bb(jind(ix),imod) = b_phi(imod,ix,ipar,jv,is)*phi(imod,ix,ipar)
          end if
          
          
          if(stress == NL_APAR) then
            aa(jind(ix),imod) = a_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
            bb(jind(ix),imod) = b_apar(imod,ix,ipar,jv,is)*apar(imod,ix,ipar)
          end if
          
          if(stress == NL_BPAR) then
            aa(jind(ix),imod) = a_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
            bb(jind(ix),imod) = b_bpar(imod,ix,ipar,jv,is)*bpar(imod,ix,ipar)
          end if
          
        end do loop_nx1
        
        ! Inverse fourier transform of field from k space to real space
        ! in radial direction (first dimension).
        ! aa = grad_zeta <field> =
        ! poloidal gradient of the gyroaverage of field in real space
        ! bb = grad_psi <field> =
        ! radial gradient of the gyroaverage of field in real space
        call fourcol(aa(:,:),FFT_INVERSE)
        call fourcol(bb(:,:),FFT_INVERSE)
        
        ! cc = grad_x f_k
        ! dd = grad_y f_k 
        loop_nvpar: do kt = 1, nvpar
        
          ! cc=grad_x_k f = radial gradient of the distribution funciton
          !                      = i J_0() k_psi fdisi
          cc(:,:) = (0.,0.)

          ! dd=grad_y_k f = zeta gradient of the distribution function
          !                     = i J_0() k_zeta fdisi
          dd(:,:) = (0.,0.)
          
          loop_nx2: do ix=1, nx
            cc(jind(ix),imod) =  ci1kxrh(ix)   * fdisi(indx(ifdis,imod,ix,ipar,jv,kt,is))
            dd(jind(ix),imod) =  ci1krho(imod) * fdisi(indx(ifdis,imod,ix,ipar,jv,kt,is))
          end do loop_nx2
        
          ! Inverse fourier transform of distribution function from k space to real space
          ! in radial direction (first dimension).
          ! cc = grad_zeta <field> =
          ! poloidal gradient of the gyroaverage of field in real space
          ! dd = grad_psi <field> =
          ! radial gradient of the gyroaverage of field in real space
          call fourcol(cc(:,:),FFT_INVERSE)
          call fourcol(dd(:,:),FFT_INVERSE)
        
          ! The zeta psi component = efun(ipar,2,1)
          ! The psi zeta component = -efun(ipar,2,1) = efun(ipar,1,2)
          ! Use the efun of the first radial grid point, since no radial
          ! dependence in flux-tube version.
          if(stress == NL_PHI) then
            dum = -efun(1,ipar,2,1)
          end if
          ! No minus sign, since apar correction is negative
          if(stress == NL_APAR) then
            dum = vpgr(ipar,jv,kt,is) * efun(1,ipar,2,1)
          end if
          if(stress == NL_BPAR) then
            dum = -efun(1,ipar,2,1)
          end if
          
          !Evaluate nonlinear term using hermitian symmetry of the spectrum
          dum_h = (0., 0.)
          loop_mrad: do j=1, mrad
            dum_h(j,imod) = dum_h(j,imod) + c1 * dum &
            & * parseval_correction(imod) * real(conjg(aa(j,imod)) * cc(j,imod) &
            & - conjg(bb(j,imod)) * dd(j,imod))
          end do loop_mrad
          
          ! Make forward fft, i.e., transform from real space in radial
          ! direction to Fourier space
          call fourcol(dum_h,FFT_FORWARD)
          
          ! nomrmalize FFT
          dum_h = dum_h / mrad
        
          !Assign the nonlinear contribution of the toroidal mode imod to iyzero-mode
          !of dt_g, which will be used later to invert Poisson equation
          loop_nx3 : do j = 1, mrad
            if (jinv(j).ne.0) then
              if (.not.(jinv(j)==ixzero.and.imod==1)) then
                dt_g(indx(ifdis,iyzero,jinv(j),ipar,jv,kt,is)) =     & 
                    &  dt_g(indx(ifdis,iyzero,jinv(j),ipar,jv,kt,is)) + dum_h(j,imod) 
              end if
            endif
          end do loop_nx3
          
        end do loop_nvpar

      end do; end do; end do
      

      ! (ii) Invert Poisson equation with dg/dt replacing distribution
      !      function in integral part.
      call calculate_fields(dt_g)
      
      
      ! (iii) Calculate the flux surface average.
      if(proc_subset(0,0,1,1,1)) then
        dum_za = (0.E0, 0.E0)
        do ix=1, nx
          do i=1, ns
          
            ! index of phi
            iih = indx(iphi,iyzero,ix,i) 
            
            ! integral in s
            dum_za(ix) = dum_za(ix) + ints(i) * dt_g(iih)
            
          end do
        end do
        
        !finish flux surface average
        call mpireduce_sum_inplace(dum_za,shape(dum_za), COMM_S_NE)
      end if
    
      ! save flux-surface-averaged part in respective term
      if(proc_subset(0,1,1,1,1)) then
        do ix=1, nx
          ctrb(ix,imod) = dum_za(ix)
        end do
      end if
    
    end do loop_nmod

  end subroutine calc_nonlin_ctrb_hermitian


  
  
  !--------------------------------------------------------------------
  !> Generic routine that outputs kx-spectra.
  !--------------------------------------------------------------------
  subroutine kx_output_array(array, lun)
  
    use grid,         only : nx
    use mpiinterface, only : root_processor
    use io,           only : xy_fmt, ascii_fmt, append_chunk
    
    complex, dimension(nx), intent(in) :: array
    integer, intent(in) :: lun

    ! do output on root_processor
    if(root_processor) then
      call append_chunk(lun, array(:), xy_fmt, ascii_fmt)
    end if
    
  end subroutine kx_output_array
  
  
  
  
  !--------------------------------------------------------------------
  !> Generic routine that outputs kxky-spectra.
  !--------------------------------------------------------------------
  subroutine kxky_output_array(array, lun)
  
    use grid,         only : nx, nmod
    use mpiinterface, only : root_processor
    use io,           only : xy_fmt, ascii_fmt, append_chunk
    
    complex, dimension(nx,nmod), intent(in) :: array
    integer, intent(in) :: lun

    ! do output on root_processor
    if(root_processor) then
      call append_chunk(lun, array, xy_fmt, ascii_fmt)
    end if

  end subroutine kxky_output_array
  
  
  

  !--------------------------------------------------------------------
  !> The routine output() should do the output to files, using the
  !> routines provided by the io module.
  !--------------------------------------------------------------------
  subroutine output()

    use mpiinterface, only : root_processor
    use control,      only : non_linear, nlapar, nlphi, nlbpar
    use control,      only : lcollisions  

    if(.not. lcalc_zonal_evo) return
    
    ! calculate zonal phi evolution
    if (lcalc_zonal_evo) then
      call calc_largestep
    end if


    ! zonal potential related output
    if(root_processor) then
    
      ! zonal potential
      call kx_output_array(zphi_kx, i_zphi_kx)
      
      ! time derivative of zonal potential
      call kx_output_array(dt_zphi_kx, i_dt_zphi_kx)
      
    end if
    
    
    ! all linear terms combinded
    call kx_output_array(zevo_lin_kx, i_zevo_lin_kx)
    
    
    ! nonlinear terms
    if(non_linear) then
          
      ! nonlinear reynolds stress
      if(nlphi) then
        call kx_output_array(zevo03_phi_kx, i_zevo03_phi_kx)
        
        if(zevo_xy_spec) then
          call kxky_output_array(zevo03_phi_kxky, i_zevo03_phi_kxky)
        end if
          
      end if
      

      ! nonlinear maxwell stress
      if(nlapar) then
      
        call kx_output_array(zevo03_apar_kx, i_zevo03_apar_kx)
        
        if(zevo_xy_spec) then
          call kxky_output_array(zevo03_apar_kxky, i_zevo03_apar_kxky)
        end if
        
      end if
      
      
      ! nonlinear correction due to parallel magnetic field compression
      if(nlbpar) then
      
        call kx_output_array(zevo03_bpar_kx, i_zevo03_bpar_kx)
        
        if(zevo_xy_spec) then
          call kxky_output_array(zevo03_bpar_kxky, i_zevo03_bpar_kxky)
        end if
        
      end if
        
    end if !non_linear
    
    
    ! individual linear terms
    if(zevo_detail) then
      ! linear term I
      call kx_output_array(zevo01_kx, i_zevo01_kx)
      
      ! linear term II
      call kx_output_array(zevo02_kx, i_zevo02_kx)
      
      ! linear term IV
      call kx_output_array(zevo04_kx, i_zevo04_kx)
      
      ! linear term V
      call kx_output_array(zevo05_kx, i_zevo05_kx)
      
      ! linear term VII (Landau damping)
      call kx_output_array(zevo07_kx, i_zevo07_kx)
      
      ! linear term VIII
      call kx_output_array(zevo08_kx, i_zevo08_kx)
      
      if(nlbpar) then
        ! linear term X
        call kx_output_array(zevo10_kx, i_zevo10_kx)
        
        ! linear term XI
        call kx_output_array(zevo11_kx, i_zevo11_kx)
      end if
    
      ! numerical dissipation
      call kx_output_array(zevo_dpar_kx, i_zevo_dpar_kx)
      call kx_output_array(zevo_dvp_kx, i_zevo_dvp_kx)
      call kx_output_array(zevo_dperp_kx, i_zevo_dperp_kx)
      
      ! collision
      if(lcollisions) then
        call kx_output_array(zevo_coll_kx, i_zevo_coll_kx)
      end if
    end if !zevo_detail
      
  end subroutine output
  
end module diagnos_zonal_evo
