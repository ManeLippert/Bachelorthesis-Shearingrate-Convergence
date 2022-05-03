!------------------------------------------------------------------------------
!> This diagnostic calculates the Reynolds and Maxwell Stress.
!>
!> \f[
!>     m_{R,s} n_{R,s} V^Ea V_Er - \tilde B_t \tilde B_r / \mu_0
!> \f]
!> where a=(t,p,r), with t the toroidal component, p the poloidal component and
!> r the radial component (actually this is only equal to what is computed in the
!> long wavelength limit and there is an explicit velocity-space integral to account
!> for gyroaveraging). The upper index on the first velocity indicates the
!> contra-variant component is used (which is what is useful for toroidal momentum):
!> note that one of the contra-variant components V^p=0 as ExB motion is perpendicular
!> to the field line.
!>
!> The normalised contravariant momentum V^t is related to the physical V^{phi}
!> via V^{phi} = (v_ref/R_ref)*( 1/(2 pi) )*V^t.
!>
!> The diagnostic performs part of an ensemble average, by spatially
!> averaging.
!>
!> *WARNING At the moment the calculations of this diagnostic
!>  are not trustworthy, because they were never used or checked.*
!>
!>  Checking this diagnostic: some direct checks have been made by initialising
!>  a single mode and directly computing the expected Reynolds stresses. Note
!>  that the trace of the Reynolds' stress in this form is V^2 ~ E^2.
!>
!>---- LIMITATIONS --------------------------------------------------
!>
!> * At the moment only components of the Reynolds Stress tensor
!>   are calculated. Maxwell stress calculation is only a stub.
!>
!------------------------------------------------------------------------------
module diagnos_stresses

  implicit none

  private

  public :: set_default_nml_values
  public :: init, bcast, check, finalize, allocate_mem
  public :: initial_output, final_output
  public :: output

  logical, save, public :: lcalc_stresses

  integer :: lun_reynolds_stress

  !Summed over species stresses
  complex :: GamR_radradS, GamR_torradS, GamR_polradS
  complex :: GamM_radradS, GamM_torradS, GamM_polradS

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    lcalc_stresses = .false.
  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast
    
    call mpibcast(lcalc_stresses,1)
    
  end subroutine bcast

  !--------------------------------------------------------------------
  !> check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()

    ! should it warn that it is not yet ready for use in global runs?

  end subroutine check

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine init(requirements)
    use mpiinterface, only : root_processor
    use diagnos_generic, only : attach_metadata_grid
    use io, only : open_real_lu, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use global, only : DISTRIBUTION!, PHI_GA_FIELD
    use diagnos_generic, only : LOCAL_DATA!, S_GHOSTCELLS
    logical, intent(inout) :: requirements(:,:)
    
    if (.not.lcalc_stresses) return

    requirements(DISTRIBUTION,LOCAL_DATA) = .true.
    ! if(some_rhostar_switch) then
    !   if (nlphi) then
    !     requirements(PHI_GA_FIELD,LOCAL_DATA) = .true.
    !     requirements(PHI_GA_FIELD,S_GHOSTCELLS) = .true.
    !   end
    ! end
    
    if(root_processor) then
      call open_real_lu('reynolds_stress', 'diagnostic/diagnos_stresses', &
         & (/ 3 /), ascii_fmt, lun_reynolds_stress)
      call attach_metadata_grid(lun_reynolds_stress, &
         & 'time', '(GamR_radradS, GamR_torradS, GamR_polradS)', ascii_fmt)
      call attach_metadata(lun_reynolds_stress, phys_unit_key, not_avail, ascii_fmt)
      call attach_metadata(lun_reynolds_stress, description_key, not_avail, ascii_fmt)
      call attach_metadata(lun_reynolds_stress, comments_key, not_avail, ascii_fmt)

    end if
  end subroutine init

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine allocate_mem()

    if (.not.lcalc_stresses) return

  end subroutine allocate_mem

  !--------------------------------------------------------------------
  !>
  !--------------------------------------------------------------------
  subroutine finalize()
    use io, only : close_lu, ascii_fmt
    use mpiinterface, only : root_processor

    if (.not.lcalc_stresses) return
    
    ! deallocate all arrays of this diagnostic
    !if(allocated(local_array_name)) deallocate(local_array_name)

    if(root_processor) then
      ! be nice and close all logical units
      call close_lu(lun_reynolds_stress, ascii_fmt)
    end if

  end subroutine finalize

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
    use io, only : xy_fmt, ascii_fmt, append_chunk
    if (.not.lcalc_stresses) return

    call calc_largestep()

    call append_chunk(lun_reynolds_stress, &
       & real((/ GamR_radradS, GamR_torradS, GamR_polradS /)), &
       & xy_fmt, ascii_fmt)

    !Output components of the Maxwell stress - To be completed.

  end subroutine output

  !--------------------------------------------------------------------------
  !> This routine calculates the Reynolds and Maxwell stress 
  !--------------------------------------------------------------------------
  subroutine calc_largestep

    use grid,           only : nx, ns, nmu, nvpar, nsp 
    use grid,           only : nmod
    use dist,           only : fdisi
    use geom,           only : ints, efun, bn, dpfdpsi, metric, signJ
    use constants,      only : pi
    use mode,           only : krho, kxrh
    use components,     only : de, mas
    use velocitygrid,   only : intmu, intvp
    use matdat,         only : get_f_from_g 
    use fields,         only : get_averaged_phi
    use mpiinterface,   only : mpiallreduce_sum_inplace
    use diagnos_generic, only : parseval_correction
    use global,         only : r_tiny
    use dist,           only : fmaxwl
    !use dist, only : iphi ! for rhostar effects

    ! integers for the loop over all grid points 
    integer :: imod, ix, i, j, k, is

    !The two tensor sums
    complex :: Erad, Epol, Etor
    complex :: Erad_cont, Epol_cont, Etor_cont

    !Different components of the Reynolds and Maxwell stresses.
    complex :: GamR_radrad, GamR_torrad, GamR_polrad
    complex :: GamM_radrad, GamM_torrad, GamM_polrad

    complex :: fdis, phi_ga
    real    :: phi2!, apa2, bpa2

    ! An array to hold the parallel derivative of the gyroaveraged potential
    ! complex :: dphi_gads(ns)

    if (.not.lcalc_stresses) return

    ! Initialize the fluxes to zero 
    GamR_radradS = 0.
    GamR_torradS = 0.
    GamR_polradS = 0.
    GamM_radradS = 0.
    GamM_torradS = 0.
    GamM_polradS = 0.

    !Calculate the stresses.
    nsp1: do is = 1, nsp

      GamR_radrad = 0.
      GamR_torrad = 0.
      GamR_polrad = 0.
      GamM_radrad = 0.
      GamM_torrad = 0.
      GamM_polrad = 0.

      ! Spatial average.
      nmod1: do imod = 1, nmod 
        nx1:  do ix = 1, nx 

          ! Integral over the velocity space 
          nmu3: do j = 1, nmu

            ! If rhostar effects are included, calculate the parallel
            ! derivative of the gyroaveraged potential
            ! if(some_rhostar_switch) then
            !   if (nlphi) call dfieldds(PHI_GA_FIELD,imod,ix,j,is,dphi_gads)
            ! end

            nvpar3: do k = 1, nvpar

              ! Do the average over the flux surface 
              ns3: do i = 1, ns

                ! the gyro-averaged potential
                phi_ga = get_averaged_phi(imod,ix,i,j,is,fdisi) 

                ! fdis is the distribution without A_par contribution  
                fdis = get_f_from_g(imod,ix,i,j,k,is,fdisi) 
                if (abs(intvp(i,j,k,is)) < r_tiny) fdis = 0. 

                !Contravariant components of velocity. Contravariant forms useful because
                !flux surface average of contravariant form is associated with momentum transport.
                !Other nice property of $\mathcal{E}^i_j$: check trace of tensor against electric field energy.
                !Contravariant form is simpler because one component (along field line) is zero whereas
                !covariant form mixes this component.
                Erad_cont = (metric(ix,i,2,1)*kxrh(ix) + metric(ix,i,2,2)*krho(imod))*signJ*pi*dpfdpsi(ix)/bn(ix,i)**2
                Etor_cont =-(metric(ix,i,1,1)*kxrh(ix) + metric(ix,i,1,2)*krho(imod))*signJ*pi*dpfdpsi(ix)/bn(ix,i)**2
                Epol_cont = 0.
                !Previously, only covariant components were calculated but flux surface 
                !averages of these are not very useful due to near-cancellation as
                !the field line is well-aligned with the toroidal direction, and this 
                !alignment varies poloidally. 
                ! The 1 in the last index of efun indicates that this
                ! is the radial component of ExB velocity
                Erad = efun(ix,i,1,1)*kxrh(ix) + efun(ix,i,2,1)*krho(imod)
                Etor = efun(ix,i,1,2)*kxrh(ix) + efun(ix,i,2,2)*krho(imod)
                Epol = efun(ix,i,1,3)*kxrh(ix) + efun(ix,i,2,3)*krho(imod)
                
                ! As far as SRG understands, the third summand comes into
                ! play if finite rhostar terms are considered:
                ! if(some_rhostar_switch) then
                !   Erad = Erad + efun(ix,i,3,1)*conjg(dphi_gads(i))*rhostar
                !   Etor = Etor + efun(ix,i,3,2)*conjg(dphi_gads(i))*rhostar
                !   Epol = Epol + efun(ix,i,3,3)*conjg(dphi_gads(i))*rhostar
                !end if

                phi2 = &
                   ! real space volume element:
                   & ints(i) &
                   ! velspace volume element (The 2*pi which is
                   ! furthermore contained in the velocity-space
                   ! Jacobian is defined into intmu) :
                   & * intvp(i,j,k,is) * intmu(j) * bn(ix,i) &
                   ! The Maxwellian background 
                   & * fmaxwl(ix,i,j,k,is) &
                   ! The binormal integral is executed using Parseval's theorem
                   & * abs(phi_ga)**2

                ! Calculate the different components of the
                ! Reynolds Stress Tensor:
                ! \Gamma_{\psi,\psi}:
                GamR_radrad = GamR_radrad + &
                   & phi2*Erad_cont*Erad*de(ix,is)*parseval_correction(imod)
                ! \Gamma_{\zeta,\psi}:
                GamR_torrad = GamR_torrad + &
                   & phi2*Etor_cont*Erad*de(ix,is)*parseval_correction(imod)
                ! \Gamma_{s,\psi}:
                ! This should be zero due to the ExB drift being perpendicular to th field line:
                ! do we want to suppress it?
                GamR_polrad = GamR_polrad + &
                     & phi2*Epol_cont*Erad*de(ix,is)*parseval_correction(imod)
                ! Next line is test code for checking trace of Reynolds stress tensor.
                !GamR_polrad = GamR_polrad + &
                !   & phi2*(Erad_cont*Erad+Etor_cont*Etor+Epol_cont*Epol)*de(ix,is)*parseval_correction(imod)
                ! Next line just test code for checking normalisation of diagnostic...
                !GamR_polrad = GamR_polrad + &
                !   & phi2*de(ix,is)*parseval_correction(imod)

                ! Calculate the different components of the
                ! Maxwell Stress Tensor:
                ! To be completed.
              end do ns3
            end do nvpar3
          end do nmu3

        end do nx1
      end do nmod1

      ! Sum the components of the Reynolds stress over the species.
      ! BFM: remove the factor of two here
      GamR_radradS = GamR_radradS + mas(is)*GamR_radrad
      GamR_torradS = GamR_torradS + mas(is)*GamR_torrad
      GamR_polradS = GamR_polradS + mas(is)*GamR_polrad

    end do nsp1

    call mpiallreduce_sum_inplace(GamR_radradS,1);
    call mpiallreduce_sum_inplace(GamR_torradS,1);
    call mpiallreduce_sum_inplace(GamR_polradS,1);
    call mpiallreduce_sum_inplace(GamM_radradS,1);
    call mpiallreduce_sum_inplace(GamM_torradS,1);
    call mpiallreduce_sum_inplace(GamM_polradS,1);

  end subroutine calc_largestep


end module diagnos_stresses
