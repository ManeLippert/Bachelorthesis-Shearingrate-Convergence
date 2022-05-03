!-------------------------------------------------------------------------------
!> The module matdat contains all the quantities associated with the matrix
!! that contains the linear terms, and sorting and compression routines.
!! 
!! The solution vector fdisi in GKW contains both the distribution function g
!! as well as all the fields. A single matrix is built for the implicit scheme, 
!! although only the distribution function is integrated forward in time. 
!! One can think of the matrix consisting of 4 sections
!! [ g^n+1 ]  = [ M  N ] [ g^n ]
!! [ 0     ]    [ P  Q ] [ phi,apar,bpar,mom ]  
!!
!! M contains the parallel derivatives, the drifts and the trapping terms 
!! N contains the source due to the Maxwell background 
!! (which is always proportional to one of the fields)
!! P is the integral part of the field equations, i.e. the integral over g 
!! Q  is the polarization (diagonal field) terms. 
!! 
!! The whole matrix [M,N,P,Q] is stored in mat, with the indices n1,n2,n3,n4 
!! delineating the different areas: M(1,n1), N(n1+1,n2), P(n2+1,n3), Q(n3+1,n4). 
!!
!! The explicit scheme does not allow mixing of any of the areas. 
!! In the non-parallel version of the code the fields are always at the end 
!! of the solution vector, although in the parallel version of the code
!! all the ghost cell buffers appear after the end of the main solution vector.
!<-------------------------------------------------------------------------------
module matdat
  use matrix_format, only : sparse_matrix
  use structures, only : matrix_element

  implicit none
  
  private

  public :: finish_matrix_section, abort_if_bad_element
  public :: put_source, add_source
  public :: get_f_from_g, connect_parallel 
  public :: pos_par_grid, connect_rad, compress_piece, get_estimated_timestep

  public :: iac, matdat_allocate, maty
  public :: mat, mat_maxwll_background, mat_poisson, mat_field_diag
  public :: matz, matg2f, source
  public :: add_element, set_indx
  public :: matm, matd, matvpd, matperpd, matcoll, matoutflow
  public :: matn, matn_e, matn_v
  public :: mat_vpar_grd_phi, mat_vd_grad_phi_fm, mat_vpar_grad_df
  public :: mat_vdgradf, mat_trapdf_4d, mat_ve_grad_fm
  public :: mat_vpar_grd_bpar, mat_vd_grad_bpar_fm

  !---------------------------------------------------------------------------
  !> The linear terms are written as a Matrix vector multiplication and the
  !> matrix is stored in the array 'mat', in compressed format. Mat is a
  !> one dimensional array with the arrays 'ii' and 'jj' determining the
  !> position of the element of 'mat' in the 2D matrix. An explicit Euler
  !> step is then given by
  !>
  !>   f_new(ii(i)) = f_new(ii(i)) + delta_time * mat(i)*f_old(jj(i))
  !>
  !> where a summation over i and j is assumed. Different parts of the
  !> matrix are determined by the integers n1,n2,n3,n4. The index function
  !> must be used to generate ii or jj for a given point in the grid.
  !---------------------------------------------------------------------------
  type(sparse_matrix) :: mat
  type(sparse_matrix) :: mat_maxwll_background

  !> matrix part P in the above picture
  type(sparse_matrix) :: mat_poisson

  type(sparse_matrix) :: mat_field_diag

  
  !> The parallel dissipation is duplicated in the matrix matd for energetics
  type(sparse_matrix) :: matd

  !> The parallel velocity dissipation term is duplicated
  !> in matrix matvpd to use it in the energetics diagnostic
  type(sparse_matrix) :: matvpd

  !> The perpendicular dissipation terms are duplicated
  !> in matrix matperpd to use it in the energetics diagnostic
  type(sparse_matrix) :: matperpd

  !> The collision term is duplicated
  !> in matrix matcoll to use it in the energetics diagnostic
  type(sparse_matrix) :: matcoll

  !> A matrix to calculate the advection at the downstream s-boundary.
  !> This is used in the energetics diagnostic to calculate outflow.
  type(sparse_matrix) :: matoutflow
  
  
  !> The vpar grad phi term is duplicated
  !> in matrix matcoll to use it in the energetics diagnostic
  type(sparse_matrix) :: mat_vpar_grd_phi
  type(sparse_matrix) :: mat_vd_grad_phi_fm
  type(sparse_matrix) :: mat_vpar_grad_df
  type(sparse_matrix) :: mat_vdgradf
  type(sparse_matrix) :: mat_trapdf_4d
  type(sparse_matrix) :: mat_ve_grad_fm
  type(sparse_matrix) :: mat_vpar_grd_bpar
  type(sparse_matrix) :: mat_vd_grad_bpar_fm
  

  !> Array matz is used for the zonal flow calculation. This array is
  !> only used for the adiabatic electron approximation when
  !> zonal_adiabatic = .true. It is always used in combination with 
  !> the matrix maty.
  !> 
  !> in the case of zonal flows the correction in the Poission equation 
  !>            Pc = n_Re exp [-E_Re] {phi} / T_Re  
  !> can be calculated in the form 
  !>            Pc_i = [ Sum_i matz_i Rho_i] / maty_i 
  !> where Rho_i represents the gyro-centre charge density. The latter 
  !> is determined by Poisson_int (for spectral_radius = .true.) 
  !>
  type(sparse_matrix), save :: matz

  !> Array that is used for the zonal flow calculation. This array is
  !> only used for the adiabatic electron approximation when
  !> zonal_adiabatic = .true.
  type(sparse_matrix), save :: maty

  !> Complex array that contains the source in the equation for the evolution
  !> of the distribution function. This source is per definition all the
  !> terms that are independent of either the distribution function or any
  !> of the fields. (neo-classical term due to the grad-B and curvature drift
  !> in the background gradients, for instance). In time integration
  !>   f_new(i) = f_new(i) + delta_time * source(i)
  complex, save, allocatable :: source(:)

  !> the matrix necessary for the
  !> transformation from the distribution g (which contains a correction due
  !> to the parallel vector potential) to the 'real' distribution f
  !> f(i) = g(i) + matg2f%mat(i)*fdis(matg2f%jj(i))
  type(sparse_matrix) :: matg2f

  !> matrix for collisions conservation 'fields'
  type(sparse_matrix) :: matm

  !> Complex array that contains the matrix elements for the calculation 
  !> of the neoclassical particle fluxes from parallel friction
  !> There is no iin because the neoclassical fluxes
  !> diagnostic is not interested in the spatial structure.
  type(sparse_matrix) :: matn
 
  !> Complex array that contains the matrix elements for the calculation 
  !> of the neoclassical heat fluxes from parallel friction
  type(sparse_matrix) :: matn_e
 
  !> Complex array that contains the matrix elements for the calculation 
  !> of the neoclassical momentum fluxes from parallel friction
  type(sparse_matrix) :: matn_v

  !> integer array for the implicit scheme
  integer, save, allocatable :: iac(:)

  integer, parameter :: ierr_UNDEFINED = -34266234
  integer, parameter :: ierr_OK         = 0
  integer, parameter :: ierr_BAD_ALL    = 11
  integer, parameter :: ierr_BAD_S      = 21
  integer, parameter :: ierr_BAD_VPAR   = 22

  !> Flag which is used to check for bad matrix elements
  logical, save :: l_bad_elem = .false.


  integer, parameter :: timestep_estimator_nmaxt = 64
  type (matrix_element), save :: timestep_estimator_max_el(0:timestep_estimator_nmaxt)
  integer, save :: timestep_estimator_nterm = 0
  
  interface compress_piece
    module procedure compress_piece_complex
    module procedure compress_piece_real
  end interface
  
contains


!-----------------------------------------------------------------------------
!> This routine allocates the arrays of the module matdat
!-----------------------------------------------------------------------------
subroutine matdat_allocate()

  use dist,    only : nsolc, ntot, nelem_nc, nelem_cc, nelem_g2f, nf
  use control, only : method, zonal_adiabatic, nlapar, neoclassics
  use grid,    only : nx, ns, nmu, nvpar, nsp, nmod 
  use general, only : gkw_abort
  use control, only : l_matd, l_matvpd, l_matperpd, l_matcoll, l_matoutflow
  use control, only : l_mat_vpar_grad_df, l_mat_vdgradf , l_mat_trapdf_4d
  use control, only : l_mat_ve_grad_fm, l_mat_vpar_grd_phi, l_mat_vd_grad_phi_fm
  use control, only : l_mat_vpar_grd_bpar, l_mat_vd_grad_bpar_fm
  use dist, only : nmat_factor_for_deriv
  use global, only : id_s, id_x, id_mod, id_vpar
  use matrix_format, only : create_matrix, matrix_format_gkwcrs

  integer :: ierr, i

  ! initialize the error code
  ierr = 0

  ! allocate the source array
  allocate(source(nsolc),stat=ierr)
  if (ierr /= 0) call gkw_abort('matdat_allocate: cannot allocate source')

  ! intialize the source to zero
  do i = 1, nsolc
    source(i) = (0.,0.)
  end do

  ! arrays used for zonal adiabatic
  if (zonal_adiabatic) then
    maty = create_matrix("diag. matrix Y for zonal adiabatic poisson calc",nx*ns,matrix_format_gkwcrs)
    matz = create_matrix("integral matrix Z for zonal adiabatic poisson calc",nx*ns,matrix_format_gkwcrs)
  end if

  mat = create_matrix("matrix of the linear terms", ntot)
  mat_maxwll_background = create_matrix("matrix of the linear Maxwell&
     & background terms", ntot)

  ! FIXME find a better estimate, this is too much
  mat_poisson = create_matrix("poisson matrix", ntot,matrix_format_gkwcrs)

  mat_field_diag = create_matrix("field matrix, diag part", ntot,matrix_format_gkwcrs)


  ! NOTE the arrays here are allocated much larger than necessary
  ! If this scheme is to stay, they should be reduced, and compression 
  ! used. 
  
  ! store the parallel dissipation in a separate matrix
  if(l_matd) then
    i = ceiling((nmat_factor_for_deriv(id_s) + 1.0)*nf)
    matd = create_matrix("parallel dissipation matrix", i, &
         & matrix_format_gkwcrs)
  end if

  ! store the parallel velocity dissipation in a separate matrix
  if(l_matvpd) then
    i = ceiling((nmat_factor_for_deriv(id_vpar) + 1.0)*nf)
    matvpd = create_matrix("parallel velocity dissipation matrix", i, &
           & matrix_format_gkwcrs)
  end if

  ! store the perpendicular dissipation in a separate matrix
  if(l_matperpd) then
    i = ceiling((nmat_factor_for_deriv(id_x) + &
       & nmat_factor_for_deriv(id_mod) + 1.0)*nf)
    matperpd = create_matrix("perpendicular dissipation matrix", i, &
             & matrix_format_gkwcrs)
  end if

  ! store the collision operator in a separate matrix
  if(l_matcoll) then
    matcoll = create_matrix("collision operator matrix", &
            & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if

  ! store the outflow terms in a separate matrix
  if(l_matoutflow) then
    i = ceiling((nmat_factor_for_deriv(id_s) + 1.0)*nf)
    matoutflow = create_matrix("outflow matrix", i,matrix_format_gkwcrs)
  end if
  
  
  ! arrays holding individual linear terms
  if(l_mat_vpar_grad_df) then
    mat_vpar_grad_df = create_matrix("I: vpar_grad_df matrix", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_vdgradf) then
    mat_vdgradf = create_matrix("II: vdgradf matrix", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_trapdf_4d) then
    mat_trapdf_4d = create_matrix("IV: trapdf_4d matrix", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_ve_grad_fm) then
    mat_ve_grad_fm = create_matrix("V: ve_grad_fm matrix", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_vpar_grd_phi) then
    mat_vpar_grd_phi = create_matrix("VII: vpar_grd_phi (Landau damping) &
      &  matrix", 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_vd_grad_phi_fm) then
    mat_vd_grad_phi_fm = create_matrix("VIII: vd_grad_phi_fm matrix", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if
  if(l_mat_vpar_grd_bpar) then
    mat_vpar_grd_bpar = create_matrix("X: Trapping due to the perturbed &
      & magnetic field", 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2), &
      & matrix_format_gkwcrs)
  end if
  if(l_mat_vd_grad_bpar_fm) then
    mat_vd_grad_bpar_fm = create_matrix("XI: vd_grad_bpar_fm", &
      & 9*nmod*nx*ns*nsp*(nmu+2)*(nvpar+2),matrix_format_gkwcrs)
  end if

  ! the arrays for the correction due to A||, i.e. for the f-to-g
  ! conversion matrix
  if (nlapar) then
    matg2f = create_matrix("f-to-g conversion matrix", nelem_g2f,matrix_format_gkwcrs)
  end if
  
  ! the arrays for collisions conservation
  if (nelem_cc > 0) then
    matm = create_matrix("collisions conservation matrix", nelem_cc,matrix_format_gkwcrs)
  end if
  

  ! the arrays for the neoclassic fluxes from parallel friction
  if (neoclassics) then
    matn = create_matrix("collisions conservation matrix", nelem_nc,matrix_format_gkwcrs)
    matn_e = create_matrix("collisions conservation matrix", nelem_nc,matrix_format_gkwcrs)
    matn_v = create_matrix("collisions conservation matrix", nelem_nc,matrix_format_gkwcrs)
  end if

  ! array used for the implicit scheme(s)
  if (method == 'IMP') then
    allocate(iac(nsolc+1),stat=ierr)
    if (ierr /= 0) then
      call gkw_abort('matdat_allocate: cannot allocate iac')
    end if
  end if  

end subroutine matdat_allocate





!-----------------------------------------------------------------------------
!> This routine stores the source
!-----------------------------------------------------------------------------
subroutine put_source(iih,mat_elem)

  integer, intent(in) :: iih
  complex, intent(in) :: mat_elem

  source(iih) = source(iih) + mat_elem

end subroutine put_source


!-----------------------------------------------------------------------------
!> Stores size of current matrix section 
!> and sorts and compress each section when complete
!> Must only be called once, sequentially, for each isel
!-----------------------------------------------------------------------------
subroutine finish_matrix_section(isel)

  use dist, only : nelem_cc
  use control, only : method, nlapar
  use control, only : l_matd, l_matvpd, l_matperpd, l_matcoll, l_matoutflow
  use control, only : l_mat_vpar_grad_df, l_mat_vdgradf , l_mat_trapdf_4d
  use control, only : l_mat_ve_grad_fm, l_mat_vpar_grd_phi, l_mat_vd_grad_phi_fm
  use control, only : l_mat_vpar_grd_bpar, l_mat_vd_grad_bpar_fm
  use general, only : gkw_abort
  use matrix_format, only : compress_matrix, finish_matrix
  integer, intent(in) :: isel

  select case(isel)

    case(3)
      call finish_matrix(mat_poisson)

      ! Check that this part of the matrix only contains field integrals
      ! if (maxval(mat%jj(n2+1:n3)) >= irs3) call gkw_abort('Bad matrix section 3')

      if (nelem_cc > 0) then
        call finish_matrix(matm)
      end if

      ! Compress the energetics matrices
      if(l_matd) then
        call finish_matrix(matd)
      end if
      if(l_matvpd) then
        call finish_matrix(matvpd)
      end if
      if(l_matperpd) then
        call finish_matrix(matperpd)
      end if
      if(l_matcoll) then
        call finish_matrix(matcoll)
      end if
      if(l_matoutflow) then
        call finish_matrix(matoutflow)
      end if
      if(l_mat_vpar_grad_df) then
        call finish_matrix(mat_vpar_grad_df)
      end if
      if(l_mat_vdgradf) then
        call finish_matrix(mat_vdgradf)
      end if
      if(l_mat_trapdf_4d) then
        call finish_matrix(mat_trapdf_4d)
      end if
      if(l_mat_ve_grad_fm) then
        call finish_matrix(mat_ve_grad_fm)
      end if
      if(l_mat_vpar_grd_phi) then
        call finish_matrix(mat_vpar_grd_phi)
      end if
      if(l_mat_vd_grad_phi_fm) then
        call finish_matrix(mat_vd_grad_phi_fm)
      end if
      if(l_mat_vpar_grd_bpar) then
        call finish_matrix(mat_vpar_grd_bpar)
      end if
      if(l_mat_vd_grad_bpar_fm) then
        call finish_matrix(mat_vd_grad_bpar_fm)
      end if

    case(4) !Final
      call finish_matrix(mat_field_diag)
      if(nlapar) then
        call finish_matrix(matg2f)
      end if

      ! Some simple tests on the matrix (could also go into add_element)
      ! if (n4 > ntot) call gkw_abort('compress_matrix: n4 too large')
      ! do i = 1, n4
      !   if ((mat%ii(i) > nsolc) .or. (mat%ii(i) < 1)) call gkw_abort('ii - err')
      !   if ((mat%jj(i) > msolc) .or. (mat%jj(i) < 1)) call gkw_abort('jj - err')
      ! end do   

      ! make checks and adjustments for implicit scheme if needed. 
      if (method =='IMP') call make_mat_imp 
            
    case default
      call gkw_abort('finish_matrix_selection: isel out of range!')
      
  end select

end subroutine finish_matrix_section


!-----------------------------------------------------------------------------
!> Compress complex matrix piece between istart and iend.
!> Same as compress_piece_real, matched via interface to compress_piece.
!>
!> iend returns the updated compressed location of the segment end.
!> Ranges must always be compressed sequentially to avoid gaps and junk data.
!> All of the matrix below istart must already have been compressed
!> It is acceptable however to recompress a range with a new piece at the end.
!>
!> If sections are compressed together after the per-section compression 
!> this results in a re-sorting but not a size reduction because there is no 
!> overlap between the four matrix sections n1,n2,n3,n4
!-----------------------------------------------------------------------------
subroutine compress_piece_complex(istart,iend,iii,jjj,mmm)

  use control, only : root_and_not_silent
  use general, only : gkw_abort

  integer, intent(in) :: istart  !< matrix position to start from
  integer, intent(inout) :: iend !< matrix position to end at, and return
  integer, intent(inout), dimension(:) :: iii !< row location data
  integer, intent(inout), dimension(:) :: jjj !< column location data
  complex, intent(inout), dimension(:) :: mmm !< matrix data

  integer :: i, ireduced, ncmp
  real :: dummatr(1:0)  !< Dummy zero sized array
  
  ! number of elements to compress          
  ncmp = iend - istart + 1

  if (ncmp < 0) call gkw_abort('Error in call to compress_piece')
  if (ncmp < 2) return ! can't sort or compress a scalar

  ! some safety checks
  if (size(iii)/=size(jjj)) call gkw_abort('compress_piece: size iii/= jjj')
  if (size(mmm)/=size(jjj)) call gkw_abort('compress_piece: size mmm/= jjj')
  if (istart > size(mmm)) call gkw_abort('compress_piece: istart fail')
  if (iend > size(mmm)) call gkw_abort('compress_piece: iend fail') 
   
  !first sort the matrix (use 6th argument to sort complex matrix)  
  !call sort_matrix(ncmp,istart,iii,jjj,dummatr,mmm) 
  !Above line has bug (see issue 126), workaround below
  call sort_matrix(ncmp,1,iii(istart:iend),jjj(istart:iend),dummatr,mmm(istart:iend)) 
     
  ireduced = istart

  ! compress the matrix (not robust if the sort was incorrect)
  do i = 1+istart, iend
    if (iii(i) == iii(ireduced) .and. jjj(i) == jjj(ireduced)) then
      mmm(ireduced) = mmm(ireduced) + mmm(i)
    else
      ireduced = ireduced + 1
      iii(ireduced) = iii(i)
      jjj(ireduced) = jjj(i)
      mmm(ireduced) = mmm(i)
    end if
  end do  

  ! Reportage
  if (root_and_not_silent) then
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*) 'Partial matrix compression successfully completed'
        write(*,223) iend,ireduced
  223  format(' Original ',I8,' elements. New ',I8,' elements')
        write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        write(*,*)
  end if
  
  iend=ireduced

end subroutine compress_piece_complex


!-----------------------------------------------------------------------------
!> Compress real matrix piece between istart and iend
!> Same as compress_piece_complex, matched via interface to compress_piece
!-----------------------------------------------------------------------------
subroutine compress_piece_real(istart,iend,iii,jjj,mmm)

  use control, only : root_and_not_silent
  use general, only : gkw_abort
  
  integer, intent(in) :: istart !< matrix position to start from
  integer, intent(inout) :: iend !< matrix posistion to end at
  integer, intent(inout), dimension(:) :: iii !< row location data
  integer, intent(inout), dimension(:) :: jjj !< column location data
  real, intent(inout), dimension(:) :: mmm !< matrix data

  integer :: i, ireduced, ncmp
  complex :: dummatc(1:0)  !< Dummy zero sized array
  
  ! number of elements to compress
  ncmp = iend - istart + 1

  if (ncmp < 0) call gkw_abort('Error in call to compress_piece')
  if (ncmp < 2) return ! can't sort or compress a scalar
  
  ! some safety checks
  if (size(iii)/=size(jjj)) call gkw_abort('compress_piece: iii/= jjj')
  if (size(mmm)/=size(jjj)) call gkw_abort('compress_piece: mmm/= jjj')
  if (istart > size(mmm)) call gkw_abort('compress_piece: istart fail')
  if (iend > size(mmm)) call gkw_abort('compress_piece: iend fail')   
    
  !first sort the matrix (use 5th argument to sort real matrix)
  !call sort_matrix(ncmp,istart,iii,jjj,mmm,dummatc) 
  !Above line has bug (see issue 126), workaround below
  call sort_matrix(ncmp,1,iii(istart:iend),jjj(istart:iend),mmm(istart:iend),dummatc) 
     
  ireduced = istart

  ! compress the matrix
  do i = 1+istart, iend
    if (iii(i) == iii(ireduced) .and. jjj(i) == jjj(ireduced)) then
      mmm(ireduced) = mmm(ireduced) + mmm(i)
    else
      ireduced = ireduced + 1
      iii(ireduced) = iii(i)
      jjj(ireduced) = jjj(i)
      mmm(ireduced) = mmm(i)
    end if
  end do

  ! Reportage
  if (root_and_not_silent) then
        write(*,*) '------------------------------------------------------'
        write(*,*) 'Partial real matrix compression successfully completed'
        write(*,229) iend,ireduced
  229  format(' Original ',I8,' elements. New ',I8,' elements')
        write(*,*) '------------------------------------------------------'
        write(*,*)
  end if
  
  iend=ireduced

end subroutine compress_piece_real


!-----------------------------------------------------------------------------
!> This routine find the elements in the matrix where each row begins, and
!> stores it in irs(:).
!> The matrix must be fully sorted and compressed.
!> The indices found are used for faster matrix-vector multiply, with
!> the explicit scheme.
!-----------------------------------------------------------------------------
subroutine find_row_index(mat,irs,rs,re)
  use general, only : gkw_abort
  !> row indices of the start and end row of the given matrix section
  integer, intent(out) :: rs, re
  !> the nth element of irs is the index of first element of the n'th
  !> matrix row with respect to the total number of elements of the
  !> matrix. The elements of the n'th matrix row are then
  !> mat(irs(n):(irs(n+1)-1)) .
  integer, dimension(:), intent(out) :: irs
  type(sparse_matrix) :: mat
  integer :: i, irow

  irs(:)=0

  if (mat%nmat == 0) then
    ! deal with the case of zero size
    rs = 1
    re = 0
  else
    re = maxval(mat%ii(1:mat%nmat))
    rs = minval(mat%ii(1:mat%nmat))  
    
    irow = rs -1
    
    ! do i = ns, ne
      ! if (mat%ii(i) /= irow) then
        ! irow = irow + 1
        ! irs(irow) = i
      ! end if
    ! end do
    
    ! index the row starts
    do i = 1, mat%nmat
        if (irow == mat%ii(i)) then
          ! Do nothing
        else if (mat%ii(i) > irow ) then
          do
            irow = irow + 1
            if (mat%ii(i) == irow) then
               irs(irow) = i
               exit
            else
               if (irow > re ) call gkw_abort('Matdat: severe irow error 1')
               irs(irow) = i ! produces an empty row
               ! this is needed in some cases (e.g. n_s_grid = 1, k=0 mode)
            end if
          end do
        else if (mat%ii(i) < irow ) then
          call gkw_abort('Matdat: severe irow error 2')
        end if
    end do    
    
    ! Extra line marking row beyond end of the matrix for convenience
    irs(irow+1) = mat%nmat+1
    
    if(maxval(irs(rs:re))>mat%nmat) then
      call gkw_abort('Matdat: Severe error in find_row_index high')
     end if
    if(minval(irs(rs:re))<1) then
      call gkw_abort('Matdat: Severe error in find_row_index low')
    end if
    
  end if

end subroutine find_row_index


!-----------------------------------------------------------------------------
!> Sort the matrix on ii and then jj. Use in-place recursive quicksort, then
!> switch to insertion sort when the partitions are of size min_size or less.
!> (min_size = 9 was found to be optimal in most cases tested)
!>
!> To use the same routine to sort both complex and real matrices requires
!> either pointers, interfaces, or redundant arguments (the choice used here) 
!> No solution is ideal
!>
!> With this choice, a dummy array argument must be passed.  
!> For a real matrix, the 5th argument is the matrix and the 6th is dummy
!> For a complex matrix, the 6th argument is the matrix and the 5th is dummy
!-----------------------------------------------------------------------------
subroutine sort_matrix(n_elem,i_start,iii,jjj,mtr,mtc)

use general, only : gkw_abort

  integer, intent(in) :: n_elem, i_start
  integer :: istart, iend
  integer, intent(inout) :: iii(:), jjj(:) !< Arrays to sort on
  logical :: lcmplx, lreal  !< switches for type of matrix to sort 
  real, dimension(:), intent(inout) :: mtr    !< real matrix to sort
  complex, dimension(:), intent(inout) :: mtc !< complex matrix to sort
    
  !> The smallest sized partition on which to perform quicksort
  integer, parameter :: min_size = 9

  ! return if nothing to sort
  if (n_elem < 2) return

  lcmplx = size(mtc)>0
  lreal = size(mtr)>0
  
  if (lcmplx.and.lreal) then
     call gkw_abort('Sort matrix called with two matrices')
  else if ( .not. (lcmplx .or.lreal)) then 
     call gkw_abort('Sort matrix called with no matrices')
  end if
  
  if (lcmplx .and. size(mtc) < n_elem) then
     call gkw_abort('sort matrix: mtc too small')
  end if
  if (lreal .and. size(mtr) < n_elem) then
    call gkw_abort('sort matrix: mtr too small')
  end if

  istart = i_start
  iend   = i_start + n_elem - 1

  call qsort(istart,iend)

  ! internal subroutines
  contains

  !---------------------------------------------------------------------------
  recursive subroutine qsort(b,e)

    integer, intent(in) :: b, e
    integer :: ind

    if (e - b > min_size) then
      ind = median_of_3_ind(b,e)
      call partition(b,e,ind)
      call qsort(b,ind-1)
      call qsort(ind+1,e)
    else
      call isort(b,e)
    end if

  end subroutine qsort
  !---------------------------------------------------------------------------
  subroutine partition(b,e,ind)

    integer, intent(in) :: b, e
    integer, intent(inout) :: ind
    integer :: ival, jval, i

    ! the pivot values
    ival = iii(ind)
    jval = jjj(ind)

    ! swap the pivot point with the end point
    call swap_elements(ind,e)

    ! put elements to the right or left of the pivot value
    ind = b
    do i = b, e-1
      if (lesseq(iii(i),jjj(i),ival,jval)) then
        call swap_elements(i,ind)
        ind = ind + 1
      end if
    end do

    ! put the pivot value back at new ind
    call swap_elements(ind,e)

  end subroutine partition
  !---------------------------------------------------------------------------
  subroutine swap_elements(i,j)

     integer, intent(in) :: i, j
     integer :: itmp
     complex :: ctmp
     real :: rtmp

     itmp = iii(i) ; iii(i) = iii(j) ; iii(j) = itmp
     itmp = jjj(i) ; jjj(i) = jjj(j) ; jjj(j) = itmp
     if (lcmplx) then 
       ctmp = mtc(i) ; mtc(i) = mtc(j) ; mtc(j) = ctmp
     end if
     if (lreal) then 
       rtmp = mtr(i) ; mtr(i) = mtr(j) ; mtr(j) = rtmp
     end if

  end subroutine swap_elements
  !---------------------------------------------------------------------------
  function median_of_3_ind(a,c)

    integer, intent(in) :: a, c
    integer :: b, median_of_3_ind

    b = (a + c) / 2
    if (larger(iii(a),jjj(a),iii(b),jjj(b))) then
      if (larger(iii(c),jjj(c),iii(a),jjj(a))) then
        median_of_3_ind = a
      else if (larger(iii(c),jjj(c),iii(b),jjj(b))) then
        median_of_3_ind = c
      else
        median_of_3_ind = b
      end if
    else if (larger(iii(b),jjj(b),iii(c),jjj(c))) then
      if (larger(iii(c),jjj(c),iii(a),jjj(a))) then
        median_of_3_ind = c
      else
        median_of_3_ind = a
      end if
    else
      median_of_3_ind = b
    end if

  end function median_of_3_ind
  !---------------------------------------------------------------------------
  subroutine isort(b,e)

    integer, intent(in) :: b, e
    integer :: ival, jval, i, j
    complex :: val = (0.0, 0.0)
    real :: rval = 0.0

    do i = b+1, e
      ival = iii(i)
      jval = jjj(i)
      if (lcmplx) val  = mtc(i)
      if (lreal) rval = mtr(i)
      ajgtval : do j = i-1, 1, -1
        if (less(ival,jval,iii(j),jjj(j))) then
          if (lcmplx) mtc(j+1) = mtc(j)
          if (lreal) mtr(j+1) = mtr(j) 
          iii(j+1)  = iii(j)
          jjj(j+1)  = jjj(j)
        else
          exit ajgtval
        end if
      end do ajgtval
      if (lcmplx) mtc(j+1) = val
      if (lreal) mtr(j+1) = rval  
      iii(j+1)  = ival
      jjj(j+1)  = jval
    end do

  end subroutine isort
  !---------------------------------------------------------------------------
  function larger(ii1,jj1,ii2,jj2)

    integer, intent(in) :: ii1, ii2, jj1, jj2
    logical :: larger

    larger = (ii1 > ii2) .or. (ii1 == ii2 .and. jj1 > jj2)

  end function larger
  !---------------------------------------------------------------------------
  function lesseq(ii1,jj1,ii2,jj2)

    integer, intent(in) :: ii1, ii2, jj1, jj2
    logical :: lesseq

    lesseq =  (ii1 < ii2) .or. (ii1 == ii2 .and. jj1 <= jj2)

  end function lesseq
  !---------------------------------------------------------------------------
  function less(ii1,jj1,ii2,jj2)

    integer, intent(in) :: ii1, ii2, jj1, jj2
    logical :: less

    less = (ii1 < ii2) .or. (ii1 == ii2 .and. jj1 < jj2)

  end function less
  !---------------------------------------------------------------------------

end subroutine sort_matrix


!------------------------------------------------------------------------------
!> This function deterimes the position of a grid point in the parallel (along
!> the magnetic field) grid. This routine also implements the open boundary
!------------------------------------------------------------------------------
function pos_par_grid(imod, ix, i, k) 

  use control,      only : spectral_radius, vp_trap, parallel_boundary_conditions
  use grid,         only : gs, n_s_grid
  use velocitygrid, only : iblow
  use mode,         only : ixplus, ixminus

  integer, intent(in)     :: imod, ix, i, k
  integer                 :: pos_par_grid 

  ! in the case of spectral_radius = .false. there are no end points 
  if (.not.spectral_radius) then 
    pos_par_grid = 0 
    return 
  endif
  
  ! For the Dirichlet boundary condition the end point schemes are not used 
  if (parallel_boundary_conditions=='Dirichlet') then 
    pos_par_grid = 0 
    return 
  endif  

  ! A trapped particle with vp_trap = 1 is never an end point 
  if (vp_trap == 1) then
    if(iblow(k) /= 0) then
      pos_par_grid = 0 
      return 
    end if
  endif 

  ! assume first that it is not an end point 
  pos_par_grid = 0

  ! use the global parameter functions to determine if there is an end point 
  if ((gs(i) .eq. 1) .and. (ixminus(imod,ix) .eq. 0) ) pos_par_grid = -2
  if ((gs(i) .eq. n_s_grid) .and. (ixplus(imod,ix) .eq. 0) ) pos_par_grid = 2
  ! the following are for 4th order boundary conditions
  if ((gs(i) .eq. 2) .and. (ixminus(imod,ix) .eq. 0) ) pos_par_grid = -1
  if ((gs(i) .eq. (n_s_grid-1)) .and. (ixplus(imod,ix) .eq. 0) ) pos_par_grid = 1

end function pos_par_grid 

!------------------------------------------------------------------------------
!> This soubroutine is used in connection with the parallel boundary
!> conditions. It determines the point in the grid to connect with
!------------------------------------------------------------------------------
subroutine connect_parallel(elem,ingrid)

  use control,      only : vp_trap
  use grid,         only : ns, nvpar, parallel_s, gs, n_s_grid
  use mode,         only : mode_box, ixplus, ixminus, iyzero
  use velocitygrid, only : iblow, ibhig
  use general,      only : gkw_abort
  use geom,         only : geom_type

  !> @param elem the matrix element 
  type (matrix_element), intent(inout) :: elem
  !> @param ingrid True if a point in the grid can be referenced
  logical, intent(out) :: ingrid

  ! Does the grid follow the trapping condition ?
  if (vp_trap .eq. 1) then
    if (parallel_s) call gkw_abort('connect_parallel: I can not deal with '//  &
        &                          'parallel_s and vp_trap = 1 yet!')
    ! is the particle trapped ?
    if (iblow(elem%k) .ne. 0) then

      ! In the grid ? 
      if ((elem%i.ge.iblow(elem%k)).and.(elem%i.le.ibhig(elem%k))) then 
        ingrid = .true.
      else 
        ingrid = .false. 
        return 
      endif

      if (elem%iloc .lt. iblow(elem%k)) then
        elem%iloc = 2*iblow(elem%k)-elem%iloc-1
        elem%kloc  = nvpar - elem%k + 1
        return
      endif
      if (elem%iloc .gt. ibhig(elem%k)) then
        elem%iloc = 2*ibhig(elem%k) - elem%iloc + 1
        elem%kloc = nvpar - elem%k + 1
        return
      endif

      ! not bouncing
      return

    endif
  endif

  ! check if the point lies on the grid
  if ((gs(elem%iloc).ge.1).and.(gs(elem%iloc).le.n_s_grid))  then
    ingrid = .true.
    return
  endif

  if (mode_box) then

    ! A 2Dimensional array of modes is used (i.e. different ix must
    ! be connected

    ! The ky = 0 mode is always periodic, i.e. (ky=0,kx=*) modes
    ! connect to themselves at parallel boundaries.

    ! In principle this
    ! should also be taken care of by ixplus / ixminus but it does not
    ! cover the case of the kx=0 mode (neoclassics only)
    if (iyzero==elem%imod .and. geom_type/='slab') then
    
      if (elem%iloc .le. 0) then
        if (parallel_s) then
          !Ghost points exist
        else
          elem%iloc = ns + elem%iloc 
        endif
        ingrid = .true.
        return
      endif
      if (elem%iloc .gt. ns) then
        if (parallel_s) then
          !Ghost points exist
        else
          elem%iloc = elem%iloc - ns
        endif
        ingrid = .true.
        return
      endif

    else

      if (elem%iloc .le. 0) then
        if (ixminus(elem%imod,elem%ix).ne.0) then
          elem%ixloc = ixminus(elem%imod,elem%ix)
          if (parallel_s) then
            !Ghost points exist
          else
            !Use of "shear-periodicity"
            elem%iloc = ns + elem%iloc
          endif
          ingrid = .true.
          return
        else
          elem%iloc   = 0
          elem%ixloc = 0
          ingrid = .false.
          return
        endif
      endif

      if (elem%iloc .gt. ns) then
        if (ixplus(elem%imod,elem%ix) .ne. 0) then
          elem%ixloc = ixplus(elem%imod,elem%ix)
          if (parallel_s) then
            !Ghost points exist
          else
            !Use of "shear-periodicity"
            elem%iloc = elem%iloc - ns
          endif
          ingrid = .true.
          return
        else
          ingrid = .false.
          return
        endif
      endif
    endif

  else !not mode_box

    ! in general the code never reaches this point, but a not-in-grid 
    ! is implemented here 
    ingrid = .false.
    return  

  endif

  ! Something went wrong if you reach this point
  call gkw_abort('Internal error in connect_parallel')

end subroutine connect_parallel


!------------------------------------------------------------------------------
!> connection routine for the radial direction. (only used if non-spectral)
!------------------------------------------------------------------------------
subroutine connect_rad(E,ingrid,bc_shift) 

  use grid,       only : n_x_grid, gx, nx, lsendrecv_x, n_s_grid
  use dist,       only : xgp => ghost_points_x
  use mode,       only : krho
  use geom,       only : alphak_xbnd, qx  
  use general,    only : gkw_abort 
  use constants,  only : ci1
  use control,    only : radial_boundary_conditions
  use components, only : rhostar 

  ! input and output variables 
  type (matrix_element), intent(inout) :: E 
  logical, intent(out)                 :: ingrid  
  real, optional, intent(in)           :: bc_shift 

  ! the ix index in the global array 
  integer :: ix_G 

  ! The end of the grid on the local processor 
  integer :: nx_loc 

  ! the additional phase shift at the boundary 
  complex :: shift 

  ! is true if the considered cell is not a boundary or is a boundary
  ! but has a nonzero value.
  ingrid = .true. 

  ! if bc_shift is present the arrays are aways considered to be of the size 
  ! n_x_grid even for parallelization in the x-direction. This is used in the 
  ! field equations 
  if (present(bc_shift)) then 
    ix_G   = E%ixloc 
    nx_loc = n_x_grid 
    shift  = bc_shift 
  else 
    ix_G   = gx(E%ixloc)
    nx_loc = nx  
    shift  = (0.,0.) 
  endif 
  
  ! The test if the grid is sufficiently large is done here 
  if (ix_G + n_x_grid .lt. 0) then 
    call gkw_abort('Radial box size too small (likely smaller than rho_i') 
  endif 
  if (ix_G - n_x_grid .gt. n_x_grid) then 
    call gkw_abort('Radial box size too small (likely smaller than rho_i') 
  endif 

  select case(radial_boundary_conditions) 
  case('periodic') 
 
    if (ix_G < 1)        E%val = E%val * exp( ci1*krho(E%imod)*(alphak_xbnd(E%iloc)-shift)) 
    if (ix_G > n_x_grid) E%val = E%val * exp(-ci1*krho(E%imod)*(alphak_xbnd(E%iloc)-shift)) 

    ! MPI communicator is periodic in X, set in grid pp(ix)%periodic=.true.
    ! with ghost cells, no need to shift the point. Therefore a remap of the points 
    ! is necessary for the case without parallelization and the field equations only 
    if ((.not.lsendrecv_x).or.present(bc_shift)) then 
      if (ix_G < 1)        E%ixloc = E%ixloc + n_x_grid 
      if (ix_G > n_x_grid) E%ixloc = E%ixloc - n_x_grid 
    endif 

  case('Dirichlet')

    if (ix_G < 1) then 
      ingrid = .false. 
      E%val   = 0
    endif 
    if (ix_G > n_x_grid) then 
      ingrid = .false. 
      E%val   = 0 
    endif 

  case('Neuslab')

    if (ix_G < 1) then 
      ingrid = .true. 
      E%ixloc = 2 - ix_G
      E%val   = E%val
    endif 
    if (ix_G > n_x_grid) then 
      ingrid = .true.
      E%ixloc = 2*nx_loc - E%ixloc 
      E%val   = E%val
    endif

  case('Neu-Dir') 

    if (ix_G < 1) then
      E%ixloc = 1-E%ixloc 
      E%iloc  = E%iloc - n_s_grid/2 
      E%val = E%val * exp(ci1*krho(E%imod)*qx(1) / (2.E0*rhostar)) 
      if (E%iloc < 1) then 
        E%iloc = n_s_grid + E%iloc      
        E%val = E%val * exp(-ci1*krho(E%imod)*qx(1) / rhostar) 
      endif 
      ingrid = .true. 
    endif 
    if (ix_G > n_x_grid) then 
      ingrid = .false. 
      E%ixloc = nx_loc 
      E%val   = 0.  
    endif 

  case default 

    call gkw_abort('Unknown radial boundary condition option')

  end select 

  ! test for the ghost cells 
  if (lsendrecv_x.and.(.not.present(bc_shift))) then 
    if (E%ixloc .gt. nx+xgp .or. E%ixloc .lt. 1-xgp) then
      write(*,*) 'ghost_points_x = ', xgp
      call gkw_abort('Not enough x ghost points for gyroaverage')
    end if  
  endif    

end subroutine connect_rad 


!-----------------------------------------------------------------------------
!> Check if the point is in the grid, then put it into the right place in the
!> relevant matrices. Return an appropriate error code.
!-----------------------------------------------------------------------------

subroutine add_element(E,ierr)

  use dist,    only : ifdis
  use general, only : gkw_abort
  use control, only : uniform_mu_grid
  use control, only : l_matd, l_matvpd, l_matperpd, l_matcoll, l_matoutflow
  use control, only : l_mat_vpar_grad_df, l_mat_vdgradf , l_mat_trapdf_4d
  use control, only : l_mat_ve_grad_fm, l_mat_vpar_grd_phi, l_mat_vd_grad_phi_fm
  use control, only : l_mat_vpar_grd_bpar, l_mat_vd_grad_bpar_fm
  use velocitygrid,   only : elem_is_on_vpar_grid
  use matrix_format,  only : put_element

  type (matrix_element), intent(in) :: E
  integer, intent(out) :: ierr

  ierr = ierr_UNDEFINED

  ! The routine below however checks if the element is in 
  ! the parallel velocity grid.
  if (E%itloc == ifdis .and. .not. elem_is_on_vpar_grid(E)) return

  if (abs(E%val) > 1.E10) then
    write(*,*) 'Very large matrix element from term ' // trim(E%term)
    l_bad_elem=.true.
  end if

  !----------------------------------------------------------------------------
  ! First test on the element of the (possible) zonal flow correction. 
  !----------------------------------------------------------------------------
  if (trim(E%term).eq.'poisson_int' .or. &
     & trim(E%term).eq.'ampere_int' .or. &
     & trim(E%term).eq.'ampere_bpar_int') then
    call put_element(mat_poisson,E)
    ierr = ierr_OK
    return
  end if

  if (trim(E%term).eq.'poisson_dia' .or. &
     & trim(E%term).eq.'poisson_bpar' .or. &
     & trim(E%term).eq.'ampere_dia') then
    call put_element(mat_field_diag,E)
    ierr = ierr_OK
    return
  end if

  if (trim(E%term).eq.'parallel velocity dissipation' .and. l_matvpd) then
    call put_element(matvpd,E)
  end if
  if (trim(E%term).eq.'hyper_disp_perp' .and. l_matperpd) then
    call put_element(matperpd,E)
  end if
  if ((uniform_mu_grid.and. trim(E%term).eq.'Collisions differential') .or.&
       &(.not.uniform_mu_grid.and. trim(E%term).eq.'Collisions differential numu')&
       ) then
    if(l_matcoll) then
      call put_element(matcoll,E)
    end if
  end if
  if (E%outflow .and. l_matoutflow) then
    call put_element(matoutflow,E)
  end if
  if (trim(E%term).eq.'parallel dissipation' .and. l_matd) then
    call put_element(matd,E)
  end if
  if (trim(E%term).eq.'VII: vpar_grd_phi (Landau damping)' .and. &
     & l_mat_vpar_grd_phi) then
    call put_element(mat_vpar_grd_phi,E)
  end if
  if(l_mat_vd_grad_phi_fm) then
    if (trim(E%term).eq.'VIII: vd_grad_phi_fm') then
      call put_element(mat_vd_grad_phi_fm,E)
    end if
  end if
  if (trim(E%term).eq.'I: vpar_grad_df' .and. l_mat_vpar_grad_df) then
    call put_element(mat_vpar_grad_df,E)
  end if
  if (trim(E%term).eq.'II: vdgradf' .and. l_mat_vdgradf) then
    call put_element(mat_vdgradf,E)
  end if
  if (trim(E%term).eq.'IV: trapdf_4d' .and. l_mat_trapdf_4d) then
    call put_element(mat_trapdf_4d,E)
  end if
  if (trim(E%term).eq.'V: ve_grad_fm' .and. l_mat_ve_grad_fm) then
    call put_element(mat_ve_grad_fm,E)
  end if
  if (trim(E%term).eq.'X: Trapping due to the perturbed magnetic &
      & field' .and. l_mat_vpar_grd_bpar) then
    call put_element(mat_vpar_grd_bpar,E)
  end if
  if (trim(E%term).eq.'XI: vd_grad_bpar_fm' .and. &
     & l_mat_vd_grad_bpar_fm) then
    call put_element(mat_vd_grad_bpar_fm,E)
  end if


  if(E%itloc == ifdis) then
    ! this is a linear term containing a distribution fluctuation
    call put_element(mat,E)
  else
    ! this is a linear term containing a field fluctuation
    call put_element(mat_maxwll_background,E)
  end if

  ! Call the timestep estimator
  call update_timestep_estimator(E)

  ierr = ierr_OK

end subroutine add_element


!-----------------------------------------------------------------------------
!> This routine stores the largest matrix elements for the timestep estimate.
!> It is called from add_element for every term put into the matrix.
!>
!> The estimate is not perfect since it only estimates von neumann stability
!> limits from the size of the matrix elements and empirical tests.
!>
!> It does not include the operator eigenvalue stability.
!>
!> To avoid searching through all the stored values repeatedly,
!> this routine first assumes that all the elements for each term are entered
!> as a contiguous block. However, this is not always true (i.e. two elements
!> are put in one loop) in which case a slower check must be made.
!-----------------------------------------------------------------------------
subroutine update_timestep_estimator(E)
  use general,        only : gkw_abort, gkw_warn
  use mpiinterface,   only : mpibcast, mpiallreduce_maxloc, mpiallreduce_max 
  use control,        only : meth, method
  use global,         only : r_tiny
  !> matrix element to store
  type (matrix_element), intent(in) :: E
  integer :: i

  ! look for the term in our list. backwards is faster here, because
  ! the new terms are appended to the list.
  do i = timestep_estimator_nterm, 1, -1
    if (trim(timestep_estimator_max_el(i)%term) == trim(E%term)) then
      if (E%ideriv /= timestep_estimator_max_el(i)%ideriv) then
        ! avoid issue 226
        call gkw_abort('%ideriv inconsistent with %term')
      end if

      if (abs(E%val) > abs(timestep_estimator_max_el(i)%val)) then
        timestep_estimator_max_el(i) = E
      end if
      return
    end if
  end do

  ! term has never been seen before.

  ! append a new term to the list
  timestep_estimator_nterm = timestep_estimator_nterm + 1
  if (timestep_estimator_nterm > timestep_estimator_nmaxt) then
    write(*,*) timestep_estimator_max_el(:)%term
    call gkw_abort('time_est: increase timestep_estimator_nmaxt and recompile')
  end if
  timestep_estimator_max_el(timestep_estimator_nterm) = E
end subroutine update_timestep_estimator


!-----------------------------------------------------------------------------
!>
!> This function returns the overall safe estimate for all terms
!> and generates a report.
!>
!-----------------------------------------------------------------------------
function get_estimated_timestep() result(time_est)
  use general, only : gkw_abort, gkw_warn
  use mpiinterface, only : mpibcast, mpiallreduce_maxloc, mpiallreduce_max 
  use mpiinterface, only : root_processor, processor_number
  use mpiinterface, only : number_of_processors, send_to_root
  use mpiinterface, only : register_tag_range, mpibarrier
  use mpicomms, only : COMM_CART
  use control, only : meth, method
  use global, only : r_tiny
  use grid, only : gs, gmu, gvpar, gx, gsp

  real :: time_est

  real :: tmax = 0., tmax1 = 0., tmax2 = 0., tmax4 = 0.

  integer :: iproc
  integer :: tag1, tag2
  integer :: i, nterm_max
  integer :: nterm_of_that_proc

  ! the range of tags to be used
  integer :: tag_range_start1, tag_range_end_inkl1
  integer :: tag_range_start2, tag_range_end_inkl2
  type(matrix_element) :: E_recvd

  call mpiallreduce_max(timestep_estimator_nterm,nterm_max,1)

  call register_tag_range(number_of_processors, &
     & tag_range_start1, tag_range_end_inkl1)
  call register_tag_range(number_of_processors*nterm_max, &
     & tag_range_start2, tag_range_end_inkl2)

  time_est = 0.

  ! reporting max term only requires a mpi operation to decide which it is
  ! iloc = maxloc(abs(timestep_estimator_max_el(1:timestep_estimator_nterm)%val)) 

  ! check for a potential problem where some processors miss a term.
  ! This should not occur if every elements adds a zero term 
  ! using the routine register_term in linear_terms.
  ! There may be a more elegant (and complicated) fix by making the 
  ! estimator match terms across processor with a lookup table

  ! first calculate the timestep for the Alven wave, or its Electro-
  ! static limit 
  call time_est_field(tmax1)

  ! skip root
  do iproc = 1, number_of_processors-1
    tag1 = tag_range_start1 + iproc-1
    if(root_processor) then
      ! agree on number of terms to compare
      call send_to_root(nterm_of_that_proc, COMM_CART, tag1)
      do i = 1, nterm_of_that_proc
        tag2 = tag_range_start2 + iproc*nterm_max + i-1
        ! recv max value of remote proc
        call send_to_root(E_recvd, COMM_CART, tag2)
        !find or append term
        call update_timestep_estimator(E_recvd)
      end do
    elseif(iproc == processor_number) then
      ! agree on number of terms to compare
      call send_to_root(timestep_estimator_nterm, COMM_CART, tag1)
      do i = 1, timestep_estimator_nterm
        tag2 = tag_range_start2 + iproc*nterm_max + i-1
        ! send local max value of the term to root
        call send_to_root(timestep_estimator_max_el(i), COMM_CART, tag2)
      end do
    end if
  end do

  ! now root should have the comprehensive list of elements with
  ! respective maximum value, for all terms on any process.
  if(root_processor) then
    do i = 1,timestep_estimator_nterm

      tmax = abs(timestep_estimator_max_el(i)%val)

      ! keep global maximum for dt_est across all terms.
      ! CFL condition in n dimensions suggests we should sum all terms
      ! for each point individually. This would require doing this 
      ! analysis after the matrix compression (possible in principle).

      ! select derivative for von neumann stability
      select case(timestep_estimator_max_el(i)%ideriv)
      case(0)    ! algebraic term, no associated timestep
        tmax = 0.0
      case(1)    ! first derivative  (1/dx in value)
        tmax1 = max(tmax,tmax1)  
      case(2,-2) ! second derivative. e.g. collisions (1/dx*dx in value) 
        tmax2 = max(tmax,tmax2)
      case(4)    ! fourth derivative, e.g. hyper-diffusion  
        tmax4 = max(tmax,tmax4)
      case default
        call gkw_warn('time_est: error in ideriv')    
      end select

      ! report timestep needed for each term (no scaling for different schemes)
      if (timestep_estimator_max_el(i)%ideriv /= 0 .and. tmax > r_tiny) then
        write(*,*) 'Von Neumann stability factor for linear term: ', &
           & trim(timestep_estimator_max_el(i)%term)
        write(*,*) 'at global location:  (imod,  ix,  i,  j,  k, isp)' 
        write(*,774) timestep_estimator_max_el(i)%imod,&
           & gx(timestep_estimator_max_el(i)%ix), &
           & gs(timestep_estimator_max_el(i)%i), &
           & gmu(timestep_estimator_max_el(i)%j), &
           & gvpar(timestep_estimator_max_el(i)%k), &
           & gsp(timestep_estimator_max_el(i)%is)
774     format('                     ',I6, I5, I4, I4, I4, I5)
        write(*,'(A,es13.5)') ' Timestep (unscaled for scheme): ', 1./tmax
        write(*,*) '------------------------------------------------------'
        write(*,*)
      end if
    end do
  end if
  call mpibcast(tmax, 1)
  call mpibcast(tmax1, 1)
  call mpibcast(tmax2, 1)
  call mpibcast(tmax4, 1)

  ! auto select time integration scheme     
  if (method=='EXP') then
    if (meth == 99) then  
      if (tmax2 > tmax1/2.4 .or. tmax4 > tmax1/2.4) then
        call gkw_warn('Diffusive auto meth 99: setting meth = -4,RKC_4')
        meth = -4
      else
        call gkw_warn('Advective auto meth 99: setting meth = 2, RK4')
        meth = 2
      end if
    else if ((tmax2 > tmax1/2.4 .or. tmax4 > tmax1/2.4) .and. meth > 0) then
      call gkw_warn('Diffusive dominated: EXP meth -4,-6,-7 should allow larger timesteps')
    end if
  end if

  ! adjust timestep for different stability boundaries RK schemes
  ! based on empirical tests of Von Neumann stability for each term
  ! independently, does not consider operator eigenvalues or combination effects
  select case(meth)
  case(1,3)  ! no knowledge here, be conservative
    tmax = max(tmax1,tmax2,tmax4)
  case(2)
    tmax = max(tmax1/2.4,tmax2,tmax4/2.4)      
  case(-4,-6,-7)
    tmax = max(tmax1,tmax2/4.0,tmax4/8.0)
  case default
    call gkw_abort('time_est: unknown meth')        
  end select

  ! set upper limit for estimated timestep (since fields not included)
  ! this avoids strange results for adiabatic electrons
  tmax=max(tmax,40.0)

  ! return the estimated timestep
  time_est = 1./tmax

end function get_estimated_timestep

!-----------------------------------------------------------------------------
!> Check across all processors if any bad elements exist, and abort if so
!> (This method could be avoided now that mpi_abort is used)
!-----------------------------------------------------------------------------
subroutine abort_if_bad_element

  use mpiinterface,  only : mpiallreduce_or
  use general,       only : gkw_abort

  logical :: labort

  call mpiallreduce_or(l_bad_elem,labort,1)

  if (labort) call gkw_abort('Bad matrix element')

end subroutine abort_if_bad_element


!-----------------------------------------------------------------------------
!> Set the row and column position of the given matrix element to the given
!> grid indices.
!> Row and column indices are set to the same values, i.e. the element is on the
!> main diagonal.
!> To make the element off-diagonal, some of the indices of the matrix element
!> must be modified after the call to this function.
!-----------------------------------------------------------------------------
subroutine set_indx(E,imod,ix,i,j,k,is)
  type (matrix_element), intent(inout) :: E
  integer, intent(in) :: imod,ix,i,j,k,is

  ! the "row index" of the element
  E%imod = imod
  E%ix   = ix
  E%i    = i
  E%j    = j
  E%k    = k
  E%is   = is

  ! here also the "column index" of the element is initialized
  E%imloc = imod
  E%ixloc = ix
  E%iloc  = i
  E%jloc  = j
  E%kloc  = k
  E%isloc = is

end subroutine set_indx



!-------------------------------------------------------------------------------
!> This routine makes the changes and checks to the matrix necessary 
!> for the implicit scheme. The matrix is resorted on jj then ii for umfpack.
!> Also  if required, the matrices for the electromagnetic response and
!> the non-adiabatic distribution are recombined with the main matrix.
!> This routine allows explicit and implicit schemes to be maintained 
!> in the same code; whilst the matrix format is setup primarily for the 
!> explicit shceme
!------------------------------------------------------------------------------- 
subroutine make_mat_imp 
  use constants, only : czero
  use control,   only : nlapar, method 
  use dist,      only : n_phi_start, nsolc, ntot, nf
  use general,   only : gkw_abort
  use global,    only : gkw_a_equal_b_accuracy
  use matrix_format, only : compress_matrix
  integer :: i, j, k, m, ierr, it, nt4, itmp 
  logical :: doeiets 
  complex, allocatable :: vec1(:), vec2(:), mth(:)
  integer, allocatable :: iiv(:), iin(:), jjn(:), iwksp(:)

  ! FIXME check this
  integer :: n4
  n4 = mat_field_diag%nmat
  
  if (method/='IMP') call gkw_abort('make_mat_imp: error') 
    
  ! allocate the test array for the implicit scheme
  allocate(iwksp(n4),stat=ierr)
  if (ierr /= 0) call gkw_abort('compress_matrix: cannot allocate iwksp')

  ! re-sort the whole matrix together 
  ! no compression should occur because the four sections are independent
  ! if compression occurs it will break loops from e.g. n2+1 to n3

  ! umfpack has the opposite definition for the
  ! matrix format. This is solved here through
  ! the exchange of ii and jj, to sort the matrix on jj first
  itmp=n4
  ! call compress_piece(1,n4,jj,ii,mat)
  call compress_matrix(mat)
  if (itmp/=n4) call gkw_abort('Element overlap error n4')
  
  ! do a test on the array
  do i = 1, nsolc
    iwksp(i) = 0
  end do
  do i = 1, n4
    iwksp(mat%ii(i)) = 1
  end do
  do i = 1, nsolc
    if (iwksp(i) == 0) then
      write (*,*) 'Element ',i,' is not referenced'
      call gkw_abort('compress_matrix: test failed')
    end if
  end do

  ! After tests the help array is no longer needed and is deallocated.
  ierr = 0
  if (allocated(iwksp)) deallocate(iwksp,stat=ierr)
  if (ierr /= 0) call gkw_abort('Cannot deallocate iwksp in compress_matrix')
    
  ! set up an ia array
  iac(1) = 1
  i = 1
  do j = 1, nsolc-1
    do while (mat%jj(i) == j) 
      i = i + 1
      if (i>n4) then
        write(*,*) 'Problem matdat'
        stop 1
      end if
    end do  
    iac(j+1) = i
  end do
  iac(nsolc+1) = n4 + 1
    
 
  if (.not.nlapar) return


  ! allocate the help vectors  !THIS IS HUGE and unecessary
  ierr = 0 
  allocate(vec1(1:nf), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate vec1 in imp_integration')
  allocate(vec2(1:nf), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate vec2 in imp_integration') 
  allocate(iiv(n_phi_start-1), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate iiv in imp_integration') 
  allocate(iin(n4), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate iin in imp_integration') 
  allocate(jjn(n4), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate jjn in imp_integration') 
  allocate(mth(n4), stat = ierr) 
  if (ierr/=0) call gkw_abort('Can not allocate mth in imp_integration') 

  nt4 = 0 

  do j = n_phi_start, nsolc 

    vec1 = 0.
    doeiets = .false.  
    it = 0
    ! The correction for the electro-magnetic field 
    if (nlapar) then 
      do i = 1, matg2f%nmat
        if (matg2f%jj(i).eq.j) then 
          vec1(i) = matg2f%mat(i)
          it = it + 1 
          iiv(it) = i 
          doeiets = .true.  
        end if 
      end do 
    end if 

    if (doeiets) then 

      vec2 = 0.E0 
      do m = 1, it 
        k = iiv(m)  
        do i = iac(k), iac(k+1)-1 
          if (mat%ii(i).lt.n_phi_start) then 
            vec2(mat%ii(i)) = vec2(mat%ii(i)) + mat%mat(i)* vec1(k)
          endif 
        end do 
      end do          
 
      do i = 1, n_phi_start-1 
        if (.not. gkw_a_equal_b_accuracy(vec2(i), czero)) then
          nt4 = nt4 + 1 
          iin(nt4) = i 
          jjn(nt4) = j 
          mth(nt4) = vec2(i) 
        endif 
      end do 

    endif  ! doeiets 

  end do 

  ! then add the new elements to the matrix mat 
  do i = 1, nt4 
    n4 = n4 + 1
    if (n4.gt.ntot) call gkw_abort('Too few elements in matdat')  
    mat%ii(n4) = iin(i)
    mat%jj(n4) = jjn(i)
    mat%mat(n4) = mth(i) 
  end do 

  ! De-allocate the help vectors 
  deallocate(vec1) 
  deallocate(vec2) 
  deallocate(iiv) 
  deallocate(iin)
  deallocate(jjn)
  deallocate(mth)

  ! and compress again ....

  ! umfpack has the opposite definition for the
  ! matrix format. This is solved here through
  ! the exchange of ii and jj, to sort the matrix on jj first
  ! call compress_piece(1,n4,jj,ii,mat)
  call compress_matrix(mat)
    
  ! set up an ia array
  iac(1) = 1
  i = 1
  do j = 1, nsolc-1
    do while (mat%jj(i) == j) 
      i = i + 1
      if (i>n4) then
        write(*,*) 'Problem matdat'
        stop 1
      end if
    end do  
    iac(j+1) = i
  end do
  iac(nsolc+1) = n4 + 1

end subroutine make_mat_imp


!------------------------------------------------------------------------------
!> This function returns the value of f (in one point) obtained from g (fdisi) 
!------------------------------------------------------------------------------
function get_f_from_g(imod,ix,i,j,k,is,gin) 
  use control,        only : nlapar
  use index_function, only : indx 
  use dist,           only : nsolc, ifdis  
  integer, intent(in) :: imod, ix, i, j, k, is 
  complex, intent(in) :: gin(nsolc)
  complex :: get_f_from_g
  integer iih 

  iih = indx(ifdis,imod,ix,i,j,k,is) 
  if (nlapar .and. matg2f%nmat > 0) then 
    get_f_from_g = gin(iih) + matg2f%mat(iih)*gin(matg2f%jj(iih))
  else 
    get_f_from_g = gin(iih) 
  endif

end function get_f_from_g


!------------------------------------------------------------------------------
!> This function returns the matrix row (ii) value of a given element 
!------------------------------------------------------------------------------
function get_elem_ii(E)

  use dist,           only : iapar, iphi, ibpar, i_mom, i_ene, ifdis 
  use dist,           only : iphi_ga, iapar_ga, ibpar_ga 
  use index_function, only : indx
  use general,        only : gkw_abort

  type (matrix_element), intent(in) :: E

  integer :: get_elem_ii

  ! Calculate iih 
  select case(E%itype)
  case(iphi,iapar,ibpar)
    get_elem_ii = indx(E%itype,E%imod,E%ix,E%i)
  case(iphi_ga,iapar_ga,ibpar_ga) 
    get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%j,E%is)
  case(i_mom,i_ene) 
    get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%is)
  case(ifdis)
    get_elem_ii = indx(E%itype,E%imod,E%ix,E%i,E%j,E%k,E%is)
  case default
    get_elem_ii = 0
    call gkw_abort('get_elem_ii: undefined itype in term '//E%term)
  end select 

end function get_elem_ii


!------------------------------------------------------------------------------
!> This function returns the matrix column (jj) value of a given element 
!------------------------------------------------------------------------------
function get_elem_jj(E)

  use dist,           only : iapar, iphi, ibpar, i_mom, i_ene, ifdis 
  use dist,           only : iphi_ga, iapar_ga, ibpar_ga 
  use index_function, only : indx
  use general,        only : gkw_abort
  
  type (matrix_element), intent(in) :: E
  
  integer :: get_elem_jj
  
  ! Calculate jjh 
  select case(E%itloc)
  case(iphi,iapar,ibpar)
    get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc)
  case(iphi_ga,iapar_ga,ibpar_ga) 
    get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%jloc,E%isloc)
  case(i_mom,i_ene) 
    get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%isloc)
  case(ifdis)
    get_elem_jj = indx(E%itloc,E%imloc,E%ixloc,E%iloc,E%jloc,E%kloc,E%isloc)
  case default
    get_elem_jj = 0
    call gkw_abort('get_elem_jj: undefined itloc in term '//E%term)
  end select 
  
end function get_elem_jj

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
subroutine add_source(rhs, deltatime)
  use dist, only : nsolc
  complex, intent(inout) :: rhs(nsolc)
  real, intent(in) :: deltatime
  integer :: i

#if defined(blas)
  call blas_axpy(nsolc, deltatime, source, 1,rhs,1)
#else
  !$omp parallel
  !$omp do private(i)
  do i = 1, nsolc
    rhs(i) = rhs(i) + deltatime*source(i)
  end do
  !$omp end do
  !$omp end parallel
#endif

end subroutine add_source

!------------------------------------------------------------------------------
!> This routine estimates the maximum timestep associated with the field 
!> and stores it in tmax = 1/ max_time_step 
!------------------------------------------------------------------------------
subroutine time_est_field(tmax)
  use components,   only : mas, de_Gx, signz, adiabatic_electrons
  use components,   only : veta_Gx
  use mpiinterface, only : root_processor, mpiallreduce_sum
  use mpiinterface, only : mpiallreduce_max
  use mpicomms,     only : COMM_SP_NE
  use grid,         only : nsp, nmod, n_x_grid, n_s_grid, parallel_sp
  use geom,         only : sgr_dist, qx, bn_G, metric_G
  use constants,    only : pi 
  use mode,         only : krho, lxinv, mode_box
  use global,       only : r_huge, r_tiny, gkw_a_equal_b_accuracy

  real, intent(out) :: tmax  

  integer :: ix, i, is 
  real    :: time_field, mer, mir, kmin2, dum 

  ! initialize 
  tmax = 0.E0 

  ! only for kinetic electrons  
  if (adiabatic_electrons) return 

  time_field = r_huge  

  do ix = 1, n_x_grid; do i = 1, n_s_grid  

    ! the mass ratios needed in the calculation
    mir = 0.E0; mer = 0.E0 
    do is = 1, nsp
      if (signz(is) < 1) then 
        mer = mas(is) / de_Gx(ix,is)
      else 
        mir = mir + mas(is)*de_Gx(ix,is) 
      endif 
    end do 

    if (parallel_sp) then 
      call mpiallreduce_sum(mir,dum,1,COMM_SP_NE)
      mir = dum 
      call mpiallreduce_max(mer,dum,1,COMM_SP_NE)
      mer = dum 
    endif 
    
    ! Estimate the minimum krho^2
    if (nmod == 1) then 
      kmin2 = krho(1)**2*metric_G(ix,i,2,2) 
    else  
      kmin2 = krho(2)**2*metric_G(ix,i,2,2)
    endif 
    ! With modebox = .false. there are cases in which 
    ! lxinv is not set
    if (mode_box) then
      if(.not. gkw_a_equal_b_accuracy(lxinv, 0.E0)) then
        kmin2 = min(2*pi*lxinv,kmin2)
      end if
    endif 
    
    ! note qx is a global array. 
    time_field = min(2.E0*pi*qx(ix)*sgr_dist*bn_G(ix,i)* &
               & sqrt(mir*(veta_Gx(ix)+kmin2*mer)),time_field) 

  end do; end do; 
 
  if (root_processor) then 
    write(*,*)
    write(*,*) 'Von Neumann stability factor for the Alfven wave (also ES limit)' 
    write(*,'(A,es13.5)') ' Timestep (unscaled for scheme): ',time_field 
    write(*,*) '------------------------------------------------------'
    write(*,*)
  endif 

  ! There are rare case for which time_field is zero. For instance for the 
  ! some neo-classical calculations that use only the (0,0) mode. In this 
  ! case the estimate is switched off 
  if (time_field < r_tiny) then 
    tmax = 0.E0 
  else 
    tmax = 1.E0/time_field 
  endif 

end subroutine time_est_field

end module matdat
