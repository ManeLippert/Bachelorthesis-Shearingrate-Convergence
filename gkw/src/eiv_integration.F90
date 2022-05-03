!-----------------------------------------------------------------------------
!> This module does no time integration but tries to find the eigenvalues
!> of the system
!> /          \    /     \ /      \
!> |\f$\lambda f\f$ |  = | A B | | f    |
!> |    0     |    | C D | | \f$ \phi \f$ |
!> \          /    \     / \      /
!> Where f denotes the distortion of the distribution function, \f$ \phi \f$
!> denotes all the relevant fields, and A, B, C, D are block matrices
!> that correspond to the operator of the gyrokinetic equation.
!> These two equations are handled here as a system
!> \f[
!> 0 = C f + D \phi
!> \f]
!> \f[
!> \lambda f = A f + B \phi
!> \f]
!> For getting the eigenvalues the libraries petsc/slepc are used.
!> Therefore we need to define the matrix in matrix free form and a
!> function that determines the action of the operator (matrix).
!> As gkw already provides this possiblity via exp_integration, we use
!> it. Otherwise it would be neccessary to invert the matrix D to solve
!> the problem
!> \f[
!> \lambda f = (A + B D^{-1} C ) f
!> \f]
!> which would be more complicated.
!>
!> These module relies on the petsc and slepc libraries as eigenvalue
!> solvers. Without these, the code will abort.
!>
!> \todo More checks to see if it works.
!>       This includes:
!>       - Check for parallel vector fields.
!>         Done a simple check, works.
!>         A check with parallel vector potential and compresssional effects
!>         did also work.
!>       - Parallelisation over x.
!>       - With spectral.
!>       - Scan over kthrho as well as rlt. Works, as long as target is not to
!>         far away.
!>       - Check if using rhs works better with EPS_CONV_ABS (EPSSetConvergenceTest(EPS eps, EPSConv conv))
!> \todo Would be better if the setting of the initial distribution would be
!>   done in init_eigenvalue.
!> \todo Simplify the call to diagnostics_naverage in eiv_diagnostics.
!> \todo Enable restarting.
!> \todo Allow setting of more than one initial vector.
!-----------------------------------------------------------------------------
module eiv_integration
#ifdef HAVE_SLEPC
#include <slepcversion.h>
#if defined(SLEPC_VERSION_LT)
#if SLEPC_VERSION_LT(3, 8, 0)
#else
#include <slepc/finclude/slepceps.h>
  USE SLEPCEPS
#endif
#endif
#endif
  use global,       only : lenswitch

  implicit none
  private

  public  :: eiv_solve, init_eiv
  public  :: eiv_integration_check_params, eiv_integration_bcast_nml
  public  :: eiv_integration_read_nml, eiv_integration_write_nml

#ifdef HAVE_SLEPC
#include <petscversion.h>
! These files have to be included, due to the slepc/petsc interface.
! \todo Check if this can be avoided by use slepceps
! Due to a interface change in 3.6, it has to be checked for the version number.
#if defined(PETSC_VERSION_LT)
#if PETSC_VERSION_LT(3, 6, 0)
#include <finclude/petscvec.h>
#else
#if PETSC_VERSION_LT(3, 8, 0)
#include <petsc/finclude/petscvec.h>
#endif
#endif
#else
#include <finclude/petscvec.h>
#endif

#include <slepcversion.h>
#if defined(SLEPC_VERSION_LT)
#if SLEPC_VERSION_LT(3, 6, 0)
#include <finclude/slepcepsdef.h>
#include <finclude/slepcsys.h>
#else
#if SLEPC_VERSION_LT(3, 8, 0)
#include <slepc/finclude/slepcepsdef.h>
#include <slepc/finclude/slepcsys.h>
#endif
#endif
#else
#include <finclude/slepcepsdef.h>
#include <finclude/slepcsys.h>
#endif


  ! The Shell Matrix, that is the dummy for the problem.
  Mat, save      :: sm
  ! The solver context.
  EPS            :: solver
  PetscInt       :: its, nconv, nest
  PetscScalar    :: eigr, eigi
  PetscReal      :: errest
  PetscScalar, save :: target_eig
  SlepcConvMonitor, save :: ctx

#else
  complex, save :: target_eig
#endif
  !> This array stores a sequence for indexing. It is used in the functions to
  !> get the (possible) eigenmode into fdisi and vice versa.
  integer, save, allocatable, dimension(:) :: sequence

  !> The maximum number of iterations, slepc will do to find the requested
  !> number of eigenvalues.
  integer, save :: max_iterations

  !> The number of eigenvalues to be computed.
  integer, save :: number_eigenvalues

  !> The number of column vectores slepc should use for the subspace.
  integer, save :: nr_column_vec

  !> Stores the routine that should be used for the matrix vector product.
  integer, save :: mat_vec_routine

    !> Stores the setting for the part of the eigenspectrum that should be searched.
  integer, save :: which_eigenvalues

  !> Counts the calls to the matrix-vector product routine.
  integer, save :: counter

  !> Tolerance given by the user, for finding an eigenmode.
  real, save    :: tolerance
  !> Tolerance for slepc for finding an eigenmode, takes trafo of eigenvalue
  !> into account.
  real, save    :: tolerance_intern

  !> Target values for growth rate and frequency.
  real, save    :: growthrate, freq


  ! Not used at the moment as, at least so far, computing a range of the spectrum
  ! only works with symmetric matrices.
  ! These two values store the upper and lower value of the part of the spectrum
  ! that should be computed.
!  real, save    :: int_lower, int_upper

  !> For checking if the input mat_vec_routine is within range.
  integer, parameter :: number_mat_vec_routines = 2

  !> Selecting which solver should be used.
  character (len = lenswitch), save :: type_solver

  !> Selecting which solver should be used.
  character (len = lenswitch), save :: type_extraction

  interface eiv_integration_write_nml
    module procedure eiv_integration_read_nml
  end interface

  ! These have been added to avoid warnings of the type "'xyz' defined
  ! but not used".
  private :: mat_vec_product_rhs, mat_vec_product_exp
  private :: comp_eigenvalues3, comp_eigenvalues4, comp_eigenvalues5
  private :: comp_eigenvalues6, comp_eigenvalues8, comp_eigenvalues9
  private :: comp_eigenvalues11, comp_eigenvalues12

contains
  !-------------------------------------------------------------------
  !> Does the necessary initializations for the module
  !> eiv_integration.
  !-------------------------------------------------------------------
  subroutine init_eiv
#ifdef HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    use control,         only : itime
    use dist,            only : nf
    use exp_integration, only : init_explicit
    use general,         only : gkw_abort
    use global,          only : root_and_verbose
    use mpiinterface,    only : root_processor

    PetscErrorCode :: ierr = 0
    PetscMPIInt    :: rank

    integer        :: i

    itime = 1
    counter = 0

    if (root_and_verbose .and. root_processor) then
      write(*,*) 'Starting initialization of the eigenvalue solver.'
    end if

    ! Initialise slepc and MPI.
    call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
    call check_error(ierr, 'SlepcInitialize')
    call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
    call check_error(ierr, 'MPI_Comm_rank')
    
    allocate(sequence(nf), stat=ierr)
    if (ierr /= 0) call gkw_abort('Eigenvalue: Could not allocate the memory &
       & for the sequence array.')
    sequence = (/ (i, i = 0 + rank*nf, (rank+1)*nf-1) /)

    ! Initialise the explicit scheme, because it is going to be used for the
    ! computation of the vector product.
    call init_explicit()

    ! For user friendliness, some conversions of the input
    ! parameters are done.

    ! This is done here, because the timestep may have been adjusted
    ! automatically just before.
    tolerance_intern = trafo_tolerance(tolerance)
    ! Convert the growth rate and freq the user has given to the
    ! values used by slepc.
    target_eig = trafo_gf_to_eiv(growthrate, freq)

    if (root_and_verbose) then
      write(*,*) 'Eigenvalue input:'
      write(*,*) cmplx(growthrate, freq), ' converted to: ', target_eig
    end if

    ! Do the initialisation of the solver (and the shell matrix).
    call init_solver(solver, sm)

#else
    use general, only : gkw_abort

    call gkw_abort('Eigenvalue integration needs the slepc and petsc' &
                      & // ' libraries.')
#endif
  end subroutine init_eiv

  !-------------------------------------------------------------------
  !> Does the necessary initializations specifically for the solver.
  !-------------------------------------------------------------------
#ifdef HAVE_SLEPC
  subroutine init_solver(eiv_solver, shell_matrix)
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    use general,         only : gkw_abort
    use global,          only : root_and_verbose
    use grid,            only : n_total
    use mpiinterface,    only : root_processor

    ! The Shell Matrix, that is the dummy for the problem.
    Mat, intent(out)          :: shell_matrix
    ! The solver context.
    EPS, intent(out)          :: eiv_solver
    PetscErrorCode :: ierr = 0
    EPSType        :: solver_type

    if (root_and_verbose .and. root_processor) then
      write(*,*) 'Starting initialization of the eigenvalue solver.'
    end if

    ! Create the shell matrix.
    call MatCreateShell(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, &
                    & n_total, n_total, 0, shell_matrix, ierr);
    call check_error(ierr, 'MatCreateShell')

    ! Set the shell matrix.
!    call MatSetFromOptions(shell_matrix, ierr)
!    call check_error(ierr, 'MatSetFromOptions')

    ! Assemble the shell matrix.
    ! This makes problems, and in petsc snes example 62 these are not used.
    ! So don't call this functions.
!    call MatAssemblyBegin(SM, ierr)
!    call check_error(ierr, MatAssemblyBegin')
!    call MatAssemblyEnd(SM, ierr)
!    call check_error(ierr, 'MatAssemblyEnd')

    call MatSetUp(shell_matrix, ierr)
    call check_error(ierr, 'MatSetUp')

    ! Set the operation of the  shell matrix.
    select case(mat_vec_routine)
    case (1)
      call MatShellSetOperation(SM, MATOP_MULT, mat_vec_product_exp, ierr);
    case (2)
      call MatShellSetOperation(SM, MATOP_MULT, mat_vec_product_rhs, ierr);
    case default
    end select
    call check_error(ierr, 'MatShellSetOperation')

    ! Create eigensolver context
    call EPSCreate(PETSC_COMM_WORLD, eiv_solver, ierr)
    call check_error(ierr, 'EPSCreate')

    ! Set operators. In this case, it is a standard eigenvalue problem.
#if SLEPC_VERSION_LT(3, 8, 0)
    call EPSSetOperators(eiv_solver, SM, PETSC_NULL_OBJECT, ierr)
#else
    call EPSSetOperators(eiv_solver, SM, PETSC_NULL_MAT, ierr)
#endif
    call check_error(ierr, 'EPSSetOperators')

    ! Set the problem type. In our case this is a non-hermitian matrix.
    call EPSSetProblemType(eiv_solver, EPS_NHEP, ierr)
    call check_error(ierr, 'EPSSetProblemType')

    call EPSSetTarget(eiv_solver, target_eig, ierr)
    call check_error(ierr, 'EPSSetTarget')

    ! Not possible (yet?) with our matrix type.
    ! call EPSSetInterval(solver, int_lower, int_upper, ierr)
    ! call check_error(ierr, 'EPSSetInterval')

    ! Tell slepc which kind of eigenvalues to search for
    select case(which_eigenvalues)
    case (1)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_LARGEST_MAGNITUDE, ierr)
    case (2)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_SMALLEST_MAGNITUDE, ierr)
    case (3)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues3, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues3, PETSC_NULL_VEC, ierr)
#endif
    case (4)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues4, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues4, PETSC_NULL_VEC, ierr)
#endif
    case (5)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues5, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues5, PETSC_NULL_VEC, ierr)
#endif
    case (6)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues6, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues6, PETSC_NULL_VEC, ierr)
#endif
    case (7)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_TARGET_MAGNITUDE, ierr)
    case (8)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues8, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues8, PETSC_NULL_VEC, ierr)
#endif
    case (9)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues9, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues9, PETSC_NULL_VEC, ierr)
#endif
    case (10)
      ! call EPSSetWhichEigenpairs(eiv_solver, EPS_ALL, ierr)
      call gkw_abort('Spectrum slicing is only available for symmetric' &
                    & // '/hermitian problems, which is here not the case.')
    case (11)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues11, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues11, PETSC_NULL_VEC, ierr)
#endif
    case (12)
      call EPSSetWhichEigenpairs(eiv_solver, EPS_WHICH_USER, ierr)
#if SLEPC_VERSION_LT(3, 8, 0)
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues12, PETSC_NULL_OBJECT, ierr)
#else
      call EPSSetEigenvalueComparison(eiv_solver, comp_eigenvalues12, PETSC_NULL_VEC, ierr)
#endif
    case default
      call gkw_abort('Unknown value for which_eigenvalues')
    end select
    call check_error(ierr, 'EPSSetEigenvalueComparison')

    ! Specify solver type
    if ('gd' == type_solver) then
      call EPSSetType(eiv_solver, EPSGD, ierr)
    else if ('jd' == type_solver) then
      call EPSSetType(eiv_solver, EPSJD, ierr)
    else if ('krylovschur' == type_solver) then
      call EPSSetType(solver, EPSKRYLOVSCHUR, ierr)
    else
      call gkw_abort('Error: Unknown type of solver: ' // type_solver // '.')
    end if
    call check_error(ierr, 'EPSSetType')

    if ('harmonic' == type_extraction) then
      call EPSSetExtraction(eiv_solver, EPS_HARMONIC, ierr)
    else if ('ritz' == type_extraction) then
      call EPSSetExtraction(eiv_solver, EPS_RITZ, ierr)
    else
      call gkw_abort('Error: Unknown type of extraction: ' &
                    & // type_extraction // '.')
    end if
    call check_error(ierr, 'EPSSetExtraction')

    ! Set tolerance(accuracy), maximum number of iterations and specify a
    ! custom routine to check for stop conditions
    call EPSSetTolerances(eiv_solver, tolerance_intern, max_iterations, ierr)
    call check_error(ierr, 'EPSSetTolerances')
#if SLEPC_VERSION_LT(3, 8, 0)
    call EPSSetStoppingTestFunction(eiv_solver, customEPSStoppingTest, &
       & PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)
#else
    call EPSSetStoppingTestFunction(eiv_solver, customEPSStoppingTest, &
       & PETSC_NULL_VEC,PETSC_NULL_FUNCTION,ierr)
#endif
    call check_error(ierr, 'EPSSetStoppingTestFunction')
    call EPSSetStoppingTest(eiv_solver, EPS_STOP_USER,ierr)
    call check_error(ierr, 'EPSSetStoppingTest')

    ! Setting the requested number of eigenvalues. The number of column
    ! vectors and the maximum dimension of projected problem should be
    ! set by slepc.
    if (0 == nr_column_vec) then
      call EPSSetDimensions(eiv_solver, number_eigenvalues, PETSC_DECIDE, &
         & PETSC_DECIDE, ierr)
    else
      call EPSSetDimensions(eiv_solver, number_eigenvalues, nr_column_vec, &
         & PETSC_DECIDE, ierr)
    end if
    call check_error(ierr, 'EPSSetDimensions')

    call SlepcConvMonitorCreate(PETSC_VIEWER_STDOUT_WORLD,            &
     &                   PETSC_VIEWER_DEFAULT,ctx,ierr)
    call check_error(ierr, 'SlepcConvMonitorCreate')
    call EPSMonitorSet(eiv_solver, EPSMONITORCONVERGED, ctx, SlepcConvMonitorDestroy, ierr)
    call check_error(ierr, 'EPSMonitorSet')

    ! Documentation pages state that this normally need not to be called manually, as EPSSolve calls it.
    call EPSSetUp(eiv_solver, ierr)
    call check_error(ierr, 'EPSSetUp')

    if (root_processor) then
      call EPSGetType(eiv_solver, solver_type, ierr)
      call check_error(ierr, 'EPSGetType')
      write(*,*) 'Using ', trim(solver_type), ' method for computing the eigenvalue.'
      write(*,*) 'Initialization of the eigenvalue solver finished.'
      write(*,*)
    end if

    call EPSView(eiv_solver, PETSC_VIEWER_STDOUT_WORLD, ierr)
    call check_error(ierr, 'EPSView')
    if (root_processor) then
      write(*,*)
    end if

  end subroutine init_solver
#endif

  !-------------------------------------------------------------------
  !> Does the cleanup for the module eiv_integration. This includes
  !> freeing the space for use objects/variables and in doing the
  !> final diagnostics.
  !-------------------------------------------------------------------
  subroutine eiv_cleanup
#ifdef HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps

    integer :: ierr

    ! Before starting the clean up, write the final output.
    ! By doing it here, it is easy to see, that exp_integration wasn't
    ! already cleaned up.
!      call diagnostic_final_output()

    if (allocated(sequence)) deallocate(sequence)

    ! Before finalizing slepc destroy the solver and the shell matrix
    call EPSDestroy(solver,ierr)
    call check_error(ierr, 'EPSDestroy')
    call MatDestroy(SM, ierr)
    call check_error(ierr, 'MatDestroy')

    call SlepcFinalize(ierr)
    call check_error(ierr, 'SlepcFinalize')
#endif
  end subroutine eiv_cleanup

  !-------------------------------------------------------------------
  !> This is the main subroutine of this module, that does the work
  !> which it was designed for: Here the eigenvalues and modes are
  !> determined.
  !> Therefore EPSSolve is used.
  !> From the solver then the requested eigenmode is got.
  !>
  !-------------------------------------------------------------------
  subroutine eiv_solve
#ifdef HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    use mpiinterface, only : root_processor
    use perform, only : perfon, perfoff
    PetscErrorCode :: ierr = 0

    call perfon('eiv_solve',2)

    call eigenvalue_main(solver)

    ! Print the solution
    ! Due to a interface change in 3.6, it has to be checked for the version
    ! number. This implicitly assumes that petsc and slepc have the same
    ! version number.
#if defined(SLEPC_VERSION_LT)
#if SLEPC_VERSION_LT(3, 6, 0)
    call EPSPrintSolution(solver, PETSC_NULL_OBJECT, ierr)
#else
#if SLEPC_VERSION_LT(3, 8, 0)
    call EPSErrorView(solver, EPS_ERROR_RELATIVE, PETSC_NULL_OBJECT, ierr)
#else
    call EPSErrorView(solver, EPS_ERROR_RELATIVE, PETSC_NULL_VIEWER, ierr)
#endif
#endif
#else
    call EPSPrintSolution(solver, PETSC_NULL_OBJECT, ierr)
#endif

    call check_error(ierr, 'EPSPrintSolution/EPSErrorView')

    call perfoff(2)

    call eiv_diagnostics()

    call eiv_cleanup()

    if (root_processor) then
      write(*,*) 'Eigenvalue solver finished'
      write(*,*)
    end if
#endif
  end subroutine eiv_solve

  !-------------------------------------------------------------------
  !> \brief Set initial value and solve an eigenvalue system.
  !>
  !> This routine will solve the given eigenvalue system.
  !> It is assumed, that the system is already initialised. Only setting of
  !> the initial eigenvector (= initial distribution) is done here.
  !>
  !> \param eiv_solver The eigenvalue solver context which has to be solved.
  !-------------------------------------------------------------------
#ifdef HAVE_SLEPC
  subroutine eigenvalue_main(eiv_solver)
    use dist,      only : fdisi, nf    ! Needed for extracting the desired eigenmode.
    use mpiinterface, only : root_processor
    use grid,      only : n_total
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    EPS,intent(inout) :: eiv_solver
    Vec               :: vector
    PetscErrorCode    :: ierr

    call VecCreateMPI(PETSC_COMM_WORLD, nf, n_total, vector, ierr)
    call check_error(ierr, 'VecCreateMPI')

    call cp_fdisi_to_vec(vector, fdisi)
    call EPSSetInitialSpace(solver, 1, vector, ierr)
    call check_error(ierr, 'EPSSetInitialSpace')

    call VecDestroy(vector, ierr)
    call check_error(ierr, 'VecDestroy')

    if (root_processor) then
      write(*,*) 'Starting computation of the eigenvalue(s).'
      write(*,*)
    end if

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     Solve the eigensystem
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call EPSSolve(eiv_solver, ierr)
    call check_error(ierr, 'EPSSolve')

  end subroutine eigenvalue_main
#endif


  !-------------------------------------------------------------------
  !> \brief Matrix multiplication routine for shell matrix.
  !>
  !> This subroutine provides the matrix multiplication algorithm for
  !> the shell matrix that is used in the modul eiv_integration.
  !> In this case the values of the input vector x are copied to the
  !> distribution function. Then a timestep is done via
  !> exp_integration::advance_large_step_explicit. These time step is the
  !> copied  to the output vector y.
  !>
  !> \param matrix The shell matrix.
  !> \param x,y Petsc vectors, the first is the input, the second the
  !>        output vector.
  !> \param ierr Petsc error code.
  !-------------------------------------------------------------------
  subroutine mat_vec_product_exp(matrix, x, y, ierr)
#ifdef HAVE_SLEPC
    use exp_integration, only : advance_large_step_explicit
    use control,         only : itime, spectral_radius
    use dist,            only : fdisi
    use fields,          only : calculate_fields, field_solve_nonspec_wrap
    use perform,         only : perfon, perfoff

    Mat, intent(in)  :: matrix
    Vec, intent(in)  :: x
    Vec, intent(inout) :: y
    PetscErrorCode :: ierr

    call perfon('mat_vec_product_exp',2)

    ! To be on the safe side: Copy the input vector to fdisi.
    call cp_vec_to_fdisi(fdisi, x)

    ! To be on the safe side: Compute the fields before making the time
    ! stepping.
    if (spectral_radius) then
      call calculate_fields(fdisi)
    else
      call field_solve_nonspec_wrap(fdisi,0.,.false.)
    endif

    ! call explicit integration, without diagnostics it small timesteps
    call advance_large_step_explicit(itime, .false.)

!    ! At least needed with the temporary MatMult, otherwise one gets
!    ! errors 'object is in the wrong state'.
!    call VecAssemblyBegin(y, ierr)
!    call VecAssemblyEnd(y, ierr)
    call cp_fdisi_to_vec(y, fdisi)

    counter = counter + 1

    call perfoff(2)

    ierr = 0

    ! Keep the compiler quiet about the unused dummy argument.
    if (.false.) write(*,*) matrix

#else
    integer        :: matrix, x, y, ierr

    ierr = 1

    ! Keep the compiler quiet.
    if (matrix > 0) continue
    if (x > 0) continue
    if (y > 0) continue
#endif
  end subroutine mat_vec_product_exp

  !-------------------------------------------------------------------------
  !> \brief Matrix multiplication routine for shell matrix.
  !>
  !> This subroutine provides the matrix multiplication algorithm for
  !> the shell matrix that is used in the modul eiv_integration.
  !> In this case the values of the input vector x are copied to the
  !> distribution function. Then a timestep is done via
  !> exp_integration::calculate_rhs. The result is then
  !> copied to the output vector y.
  !>
  !> \param matrix The shell matrix.
  !> \param x,y Petsc vectors, the first is the input, the second the
  !>        output vector.
  !> \param ierr Petsc error code.
  !-------------------------------------------------------------------------
  subroutine mat_vec_product_rhs(matrix, x, y, ierr)
#ifdef HAVE_SLEPC
    use control,         only : spectral_radius, time
    use dist,            only : fdisi, rhsk
    use exp_integration, only : calculate_rhs
    use fields,          only : calculate_fields, field_solve_nonspec_wrap
    use fields,          only : calculate_cons
    use perform,         only : perfon, perfoff
    use rotation,        only : shear_shift_ky, shear_ky_shift
    use rotation,        only : shear_remap, wavevector_remap

    Mat, intent(in)  :: matrix
    Vec, intent(in)  :: x
    Vec, intent(inout) :: y
    PetscErrorCode   :: ierr

    call perfon('mat_vec_product_rhs',2)

    ! To be on the safe side: Copy the input vector to fdisi.
    call cp_vec_to_fdisi(fdisi, x)

    ! To be on the safe side: Compute the fields before making the time
    ! stepping.
    if (spectral_radius) then
      call calculate_fields(fdisi)
    else
      call field_solve_nonspec_wrap(fdisi,0.,.false.)
    endif

    if (shear_ky_shift) call shear_shift_ky(fdisi)
    if (shear_remap)    call wavevector_remap(fdisi,time)

    call calculate_rhs(fdisi, rhsk(:,1))

    call calculate_cons(rhsk(:,1))

!    ! At least needed with the temporary MatMult, otherwise one gets
!    ! errors 'object is in the wrong state'.
!    call VecAssemblyBegin(y, ierr)
!    call VecAssemblyEnd(y, ierr)
    call cp_fdisi_to_vec(y, rhsk(:,1))

    counter = counter + 1

    call perfoff(2)

    ierr = 0

    ! Keep the compiler quiet about the unused dummy argument.
    if (.false.) write(*,*) matrix
#else
    integer, intent(in) :: matrix, x, y
    integer, intent(inout) :: ierr

    ierr = 1

    if (.false.) write(*,*) matrix, x, y, ierr
#endif
  end subroutine mat_vec_product_rhs

  !-------------------------------------------------------------------
  !> \brief Copy values of petsc vector to an array.
  !>
  !> This subroutine copies values from a petsc vector to an array.
  !> \note It is not checked if the size of the array is big enough
  !>       to store all the data from the vector.
  !>
  !> \param fdisi The array to which to copy the values of the vector.
  !> \param vector Petsc vector from which to get the values that
  !>        should be copied.
  !-------------------------------------------------------------------
  subroutine cp_vec_to_fdisi(fdisi, vector)
#ifdef HAVE_SLEPC
    Vec, intent(in)       :: vector
    PetscErrorCode        :: ierr
    complex, intent(out)  :: fdisi(:)
    integer               :: istart, iend

    call VecGetOwnershipRange(vector, istart, iend, ierr)
    call check_error(ierr, 'VecGetOwnershipRange')

    call VecGetValues(vector, iend-istart, sequence, fdisi, ierr)
    call check_error(ierr, 'VecGetValues')
#else
    integer :: vector, fdisi

    ! Keep the compiler quiet.
    if(vector > 0) continue
    if(fdisi > 0) continue
#endif
  end subroutine cp_vec_to_fdisi

  !-------------------------------------------------------------------
  !> \brief Copy values of petsc vector to an array.
  !>
  !> This subroutine copies values from an array to a petsc vector.
  !> \note It is not checked if the size of the vector is large enough
  !>       to store all the data from the array.
  !>
  !> \param vector Petsc vector to which to copy the values.
  !> \param fdisi The array from which to get the values.
  !-------------------------------------------------------------------
  subroutine cp_fdisi_to_vec(vector, fdisi)
#ifdef HAVE_SLEPC
    !> \todo Make subroutine work for multiple threads.

    Vec, intent(inout)    :: vector
    PetscErrorCode        :: ierr
    complex, intent(in)   :: fdisi(:)
    integer               :: istart, iend

    call VecGetOwnershipRange(vector, istart, iend, ierr)
    call check_error(ierr, 'cp_fdisi_to_vec::VecGetOwnershipRange')

    call VecSetValues(vector, iend-istart, sequence, fdisi, INSERT_VALUES, ierr)
    call check_error(ierr, 'cp_fdisi_to_vec::VecSetValues')

    ! At least needed with the temporary MatMult, otherwise one gets
    ! errors 'object is in the wrong state'.
    call VecAssemblyBegin(vector, ierr)
    call check_error(ierr, 'cp_fdisi_to_vec::VecAssemblyBegin')
    call VecAssemblyEnd(vector, ierr)
    call check_error(ierr, 'cp_fdisi_to_vec::VecAssemblyEnd')
#else
    integer :: vector, fdisi

    ! Keep the compiler quiet.
    if(vector > 0) continue
    if(fdisi > 0) continue
#endif
  end subroutine cp_fdisi_to_vec

  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues3(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,growthrate_a > growthrate_b)

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues3

  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues5(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,freq_a > freq_b)

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues5

  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues4(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,growthrate_a < growthrate_b)

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues4

  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues6(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,freq_a < freq_b)

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues6

  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues8(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,abs(growthrate_a - growthrate) < abs(growthrate_b - growthrate))

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues8


  !-------------------------------------------------------------------
  !>
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues9(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    ! ternary operator
    res = merge(-1,1,abs(freq_a - freq) < abs(freq_b - freq))

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues9

  !-------------------------------------------------------------------
  !> Subroutine to be called if the programmer defined comparison of
  !> the eigenvalues should be used.
  !>
  !> This comparison routine simply prefers the eigenvalue with the
  !> smaller positive real part. If the real parts are the same, or
  !> both are negative, none is preferred.
  !>
  !> \note As petsc should be configured with using complex number,
  !>       the parameters ai and bi are zero as the complete complex
  !>       number could be stored in ar/br respectively.
  !>
  !> \param ar real part of the 1st eigenvalue
  !> \param br real part of the 2nd eigenvalue
  !> \param ai, bi imaginary part of the 1st,2nd eigenvalue. They are not
  !>        used, because scalars are generally complex.
  !> \param res output parameter for the result of the comparison.
  !>        Due to the interface it has to be less than zero if the
  !>        first eigenvalue is preferred, greater than zero if the
  !>        second eigenvalue is preferred and zero if none is more
  !>        preferred than the other.
  !> \param ctx input integer for additional information, optional context
  !> \param ierr output parameter for the error state.
  !>
  !> The returned parameter res can be
  !>    negative  - if the 1st eigenvalue is preferred to the 2st one
  !>    zero      - if both eigenvalues are equally preferred
  !>    positive  - if the 2st eigenvalue is preferred to the 1st one
  !>
  !> \note The interface of this subroutine is determined by SLEPc and must
  !>   not be changed.
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues11(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr

    ! Real variables for the comparison.
    real :: growthrate_a, growthrate_b
    real :: freq_a, freq_b

    ierr = 0

    call trafo_eiv_to_gf(growthrate_a, freq_a, ar)
    call trafo_eiv_to_gf(growthrate_b, freq_b, br)

    if((growthrate_a >= 0) .eqv. (growthrate_b >= 0)) then
      ! if they have the same sign
      res = (abs(growthrate_a) - abs(growthrate_b))/abs(abs(growthrate_a) - abs(growthrate_b))
    else
      ! ternary operator
      res = merge(-1,1,growthrate_a > growthrate_b)
    end if

    ! To silence compiler warnings about unused parameters.
    if (.false.) res = int(ar * ai * br * bi * ctx)
  end subroutine comp_eigenvalues11

  !-------------------------------------------------------------------
  !> Subroutine to be called if the programmer defined comparison of
  !> the eigenvalues should be used.
  !>
  !> This comparison routine prefers the eigenvalue which is closer
  !> to the target value in the complex plane.
  !>
  !> \note As petsc should be configured with using complex number,
  !>       the parameters ai and bi are zero as the complete complex
  !>       number could be stored in ar/br respectively.
  !>
  !> \param ar first eigenvalue.
  !> \param br second eigenvalue.
  !> \param ai, bi are not used (because of using complex datatype).
  !> \param res output parameter for the result of the comparison.
  !>        Due to the interface it has to be less than zero if the
  !>        first eigenvalue is preferred, greater than zero if the
  !>        second eigenvalue is preferred and zero if none is more
  !>        preferred than the other.
  !> \param ctx input integer for additional information.
  !> \param ierr output parameter for the error state. 
  !-------------------------------------------------------------------
  subroutine comp_eigenvalues12(ar, ai, br, bi, res, ctx, ierr)
    complex, intent(in)  :: ar, ai, br, bi
    integer, intent(out) :: res
    integer, intent(in)  :: ctx
    integer, intent(out) :: ierr
    ! Real variables for the comparison :
    ! Distances of a/b to the target point.
    real                 :: absa, absb

    ierr = 0

    absa = abs(ar - target_eig)
    absb = abs(br - target_eig)

    if (absa > absb) then
      res =  1
    else if (absa < absb) then
      res = -1
    else
      res =  0
    end if

    ! To silence compiler warnings about unused parameters.
    ! False is used instead of logical_false to not add dependencies to this
    ! subtroutine.
    !> \note The interface of this subroutine is determined by SLEPc and must
    !>   not be changed.
    if (.false.) res = int(ar * ai * br * bi * ctx)

  end subroutine comp_eigenvalues12

  !-------------------------------------------------------------------
  !> Subroutine for reading and writing the namelist of the module.
  !>
  !> \param ilun device number of the file for reading/writing the
  !>             namelist.
  !> \param io_stat return value for success/error.
  !> \param lwrite  optional logical, if true or not present, the name-
  !>                list is written instead of read from the file.
  !-------------------------------------------------------------------
  subroutine eiv_integration_read_nml(file_unit, io_stat, lwrite)
#ifdef HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps, only : PETSC_DECIDE
#endif
    use mpiinterface, only : root_processor
    use io, only : write_run_parameter
    integer, intent(in)  :: file_unit
    integer, intent(out) :: io_stat
    logical, optional, intent(in) :: lwrite

    ! deprecated, do not use these. They are only initialised to make
    ! valgrind quiet.
    integer :: comparison_routine
    logical :: luse_initial_value

    namelist /eiv_integration/ max_iterations, tolerance, number_eigenvalues, &
             mat_vec_routine, comparison_routine, which_eigenvalues, &
             growthrate, freq, luse_initial_value, type_solver, type_extraction, &
             nr_column_vec

    luse_initial_value = .true.
    comparison_routine = 1
    io_stat = 0

    if (present(lwrite)) then

      if (.not. lwrite) then

        ! Before reading the namelist set the default values, to be
        ! sure the variable is initialized at the end of the routine.
#ifdef HAVE_SLEPC
        ! This parameter does not correspond to a sensible number
        max_iterations = PETSC_DECIDE
#else
        max_iterations = 1
#endif
        number_eigenvalues = 1
        ! This is later handled as PETSC_DECIDE in the initialization routine:
        nr_column_vec      = 0
        mat_vec_routine    = 1   ! Using advance_large_step_explicit
        which_eigenvalues  = 11  ! user defined.

        tolerance          = 1.0e-4
        growthrate         = 1.0e+0
        freq               = 0.0e+0
!        int_lower          = 0.001
!        int_upper          = 0.500

        type_solver        = 'krylovschur'
        type_extraction    = 'ritz'

        read(file_unit, NML=eiv_integration, IOSTAT=io_stat)

      end if

    else
      ! write to input.out if called without the switch; this is the default

      ! write the namelist
      if(root_processor) write(file_unit, NML=eiv_integration)

      call write_run_parameter('eiv_integration', 'max_iterations', max_iterations)
      call write_run_parameter('eiv_integration', 'tolerance', tolerance)
      call write_run_parameter('eiv_integration', 'number_eigenvalues', number_eigenvalues)
      call write_run_parameter('eiv_integration', 'mat_vec_routine', mat_vec_routine)
      call write_run_parameter('eiv_integration', 'which_eigenvalues', which_eigenvalues)
      call write_run_parameter('eiv_integration', 'growthrate', growthrate)
      call write_run_parameter('eiv_integration', 'freq', freq)
      call write_run_parameter('eiv_integration', 'type_solver', type_solver)
      call write_run_parameter('eiv_integration', 'type_extraction', type_extraction)
      call write_run_parameter('eiv_integration', 'nr_column_vec', nr_column_vec)
    end if
  end subroutine eiv_integration_read_nml

  !-------------------------------------------------------------------
  !> check if the settings for this module are ok and aborts if not.
  !-------------------------------------------------------------------
  subroutine eiv_integration_check_params
#ifdef HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps, only : PETSC_DECIDE
#endif
    use control,      only : method, normalized
    use general,      only : gkw_abort, gkw_warn

    logical           :: abort = .false.

    ! If the eigenvalue solver should not be used, there is no need to
    ! to any checks.
    ! This might make problems if another module requires the eigenvalue
    ! solver.
    if(method /= 'EIV') return

    if ((mat_vec_routine < 1) .or. &
      & (number_mat_vec_routines > mat_vec_routine)) then
    end if

    if ((which_eigenvalues < 1) .or. (which_eigenvalues > 12)) then
      call gkw_warn('eiv_integration: which_eigenvalues out of range (1-12).')
      abort = .true.
    else if (which_eigenvalues == 10) then
      call gkw_abort('eiv_integration: which_eigenvalues can not have value 10'&
        &//'(EPS_ALL), as this is not applicable to this problem.')
      abort = .true.
    end if

    if (tolerance <= 0) then
      call gkw_warn('eiv_integration: tolerance has to be greater than zero.')
      abort = .true.
    end if

    if (max_iterations <= 0 &
#ifdef HAVE_SLEPC
       & .and. max_iterations /= PETSC_DECIDE &
#endif
      &) then
      call gkw_warn('eiv_integration: max_iterations has to be greater than zero.')
      abort = .true.
    end if

    if (0 >= number_eigenvalues) then
      call gkw_warn('eiv_integration: number_eigenvalues has to be greater than zero.')
      abort = .true.
    end if

    if (0 > nr_column_vec) then
      call gkw_warn('eiv_integration: nr_column_vec has to be greater than or equal to zero.')
      abort = .true.
    end if

    if(.not.normalized) then
      call gkw_warn('Please note that the value of CONTROL.normalize will be &
         &ignored and T is forced.')
      ! Normalisation is useful because the modes found have arbitrary
      ! amplitude and phase.
      normalized = .true.
    end if

#ifdef HAVE_SLEPC
#else
    call gkw_warn('Eigenvalue solver needs the petsc/slepc libraries for working.')
!    abort = .true.
#endif

    if (abort) then
      call gkw_abort('eig_integration: Found errors while checking params &
        &(see above), aborting.')
    end if
  end subroutine eiv_integration_check_params

  !-------------------------------------------------------------------
  !> This subroutine makes the read setting of this module available
  !> for all mpi threads.
  !-------------------------------------------------------------------
  subroutine eiv_integration_bcast_nml

    use mpiinterface, only : mpibcast

    call mpibcast(max_iterations,     1)
    call mpibcast(number_eigenvalues, 1)
    call mpibcast(nr_column_vec,      1)
    call mpibcast(which_eigenvalues,  1)

    call mpibcast(mat_vec_routine,    1)
    call mpibcast(tolerance,          1)
    call mpibcast(tolerance_intern,   1)
    call mpibcast(target_eig,         1)
    call mpibcast(growthrate,         1)
    call mpibcast(freq,               1)
!    call mpibcast(int_lower,          1)
!    call mpibcast(int_upper,          1)

    call mpibcast(type_solver,     lenswitch)
    call mpibcast(type_extraction, lenswitch)

  end subroutine eiv_integration_bcast_nml

  !-------------------------------------------------------------------
  !> \brief Wrapper for doing the diagnostics with the eigenvalue solver.
  !>
  !> This subroutine is a wrapper, that allows working of the diagnostics
  !> also for the eigenvalue solver.
  !> Basically the number of the converged eigenvalue is treated as
  !> timestep, for doing the diagnostics.
  !-------------------------------------------------------------------
  subroutine eiv_diagnostics
#ifdef  HAVE_SLEPC
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    use control,         only : itime, dtim, naverage
    use control,         only : spectral_radius, stop_filename, io_legacy
    use control,         only : time, last_largestep_time, last_smallstep_time
    use diagnostic,      only : diagnostic_naverage, diagnostic_final_output
    use diagnos_generic, only : lwrite_output1
    use dist,            only : fdisi, nsolc, nf
    use exp_integration, only : advance_large_step_explicit
    use fields,          only : calculate_fields, field_solve_nonspec_wrap
    use grid,            only : n_total
    use io,              only : open_real_lu, close_lu, append_chunk
    use io,              only : ascii_fmt
    use normalise,       only : normalise_fdisi, rotate_fdis
    use normalise, only : normalise_after_timestep
    use mpiinterface,    only : root_processor
    use restart,         only : write_restart_file
    use global, only : dotdat

#if defined(PETSC_USE_FORTRAN_DATATYPES)
    type(Vec)      :: vreal
    type(Vec)      :: vimag
    type(PetscErrorCode) :: ierr
#else
    Vec            :: vreal
    Vec            :: vimag
    PetscErrorCode :: ierr
#endif

    integer        :: ii, number_converged
    integer        :: reason_converged, iterations
    complex        :: eigr, eigi
    real           :: g, f, error, errest, residual
    integer        :: i_eiv

    call EPSGetConvergedReason(solver, reason_converged, ierr)

    if (root_processor) then
      ! Just printing the number wouldn't tell the user much.
      select case(reason_converged)
      case (EPS_CONVERGED_TOL)
        write(*,*) 'Solver stopped because solution converged within requested&
           & tolerance.'
      case (EPS_DIVERGED_ITS)
        write(*,*) 'Solver aborted due to reaching maximum number of iterations.'
      case (EPS_CONVERGED_USER)
        write(*,*) 'Solver terminated due to reaching the &
           &maximum runtime (max_seconds) or because '//stop_filename//' was present.'
      case (EPS_DIVERGED_BREAKDOWN)
        write(*,*) 'Solver aborted due to method breakdown.'
      case default
        write(*,*) 'Unknown reason for stopping of the solver.'
      end select
    end if

    ! Create the vectors for storing the eigenvalues.
    call VecCreateMPI(PETSC_COMM_WORLD, nf, n_total, vreal, ierr)
    call check_error(ierr, 'VecCreateMPI::vreal')
    call VecCreateMPI(PETSC_COMM_WORLD, nf, n_total, vimag, ierr)
    call check_error(ierr, 'VecCreateMPI::vimag')

    call EPSGetConverged(solver, number_converged, ierr)
    call check_error(ierr, 'EPSGetConverged')
    if (root_processor) then
      write(*,*) 'Number of converged eigenmodes: ', number_converged
    end if

    if (root_processor) then
      call EPSGetIterationNumber(solver, iterations, ierr);
      call check_error(ierr, 'EPSGetIterationNumber')
      write(*,*) 'Iterations needed: ', iterations
      write(*,*) 'Calls to mat_vec_product:', counter
      write(*,*)
    end if

    if(root_processor .and. (0 < number_converged)) then
      call open_real_lu(dotdat('eigenvalues',io_legacy), '/diagnostic', (/ 3 /), &
             & ascii_fmt, i_eiv)
    end if

    do ii = 0, number_converged - 1
      call EPSGetEigenvalue(solver, ii, eigr, eigi, ierr)
      call check_error(ierr, 'EPSGetEigenvalue')
      residual = -1.0
#if defined(SLEPC_VERSION_LT)
#if SLEPC_VERSION_LT(3, 6, 0)
      call EPSComputeRelativeError(solver, ii, error, ierr)
      call check_error(ierr, 'EPSComputeRelativeError')
      call EPSComputeResidualNorm(solver, ii, residual, ierr)
      call check_error(ierr, 'EPSComputeResidualNorm')
#else
      ! Other possibilities for the error type are EPS_ERROR_ABSOLUTE
      ! and EPS_ERROR_BACKWARD.
      call EPSComputeError(solver, ii, EPS_ERROR_RELATIVE, error, ierr)
      call check_error(ierr, 'EPSComputeError')
#endif
#else
      call EPSComputeRelativeError(solver, ii, error, ierr)
      call check_error(ierr, 'EPSComputeRelativeError')
      call EPSComputeResidualNorm(solver, ii, residual, ierr)
      call check_error(ierr, 'EPSComputeResidualNorm')
#endif
      call EPSGetErrorEstimate(solver, ii, errest, ierr)
      call check_error(ierr, 'EPSGetErrorEstimate')

      if(root_processor) then
        write(*,'(A,i5)') 'eigenmode number: ', ii+1

        call trafo_eiv_to_gf(g, f, eigr)
        write(*,*) 'slepc eigenvalue: ', eigr
        write(*,*) '->  growth rate ', g, ' freq.', f
        if(residual >= 0.0) then
          write(*,*) ' with error (est) = ', error, ' (', errest,  &
             & ') residual = ', residual
        else
          write(*,*) ' with error (est) = ', error, ' (', errest, ')'
        end if
        write(*,*)
        call append_chunk(i_eiv, &
           & (/ real(ii+1), g, f /), &
           & '(3(es13.5))', ascii_fmt)
      end if

      call EPSGetEigenvector(solver, ii, vreal, vimag, ierr)
      call check_error(ierr, 'EPSGetEigenvector')
      call cp_vec_to_fdisi(fdisi, vreal)

      ! For the eigenmode obtained from the solver, distribution
      ! and potential may not be consistent, so do this again.
      if (spectral_radius) then
        call calculate_fields(fdisi)
      else
        call field_solve_nonspec_wrap(fdisi,0.,.false.)
      endif

      call normalise_fdisi(fdisi, nsolc)
      ! Up to here, the phase of the eigenmode is arbitrary and may
      ! depend e.g. on the parallelisation. Therefore rotate in the
      ! complex plane to output a better defined eigenmode.
      call rotate_fdis(fdisi)


      ! In the following lines, the explicit solver is called a few
      ! times, just to be able to compute all relevant diagnostics
      ! (some still insist on calculating time derivatives even for
      ! converged linear modes...)

      ! These calls are for getting the correct values growth rate
      ! and frequency in time.dat
      time = real(ii+1) - 2*naverage*dtim
      last_largestep_time = time
      last_smallstep_time = time
      call advance_large_step_explicit(itime, .true.)
      call normalise_after_timestep()
      ! We are not yet interested in timedependent diagnostics, thus
      ! disable output
      lwrite_output1 = .false.
      call diagnostic_naverage()
      last_largestep_time = time
      call advance_large_step_explicit(itime, .true.)
      call normalise_after_timestep()

      ! To write the desired values to file
      lwrite_output1 = .true.
      ! After 2 large timesteps, we have now time = real(ii+1),
      ! i.e. the mode number instead of a sensible time. This is then
      ! output via the dominant growth rate diagnostic.
      call diagnostic_naverage()
      ! pass the eigenmode number through the argument,
      ! which is used to number the outputs.
      call diagnostic_final_output(ii+1)

      ! Do the final diagnostics for every
      ! mode

      call write_restart_file(.true., ii+1,.false.)
    end do

    call VecDestroy(vreal, ierr)
    call check_error(ierr, 'eiv_diagnostics::VecDestroy::vreal')
    call VecDestroy(vimag, ierr)
    call check_error(ierr, 'eiv_diagnostics::VecDestroy::vimag')

    if(root_processor .and. (0 < number_converged)) then
      call close_lu(i_eiv, ascii_fmt)
    end if

#endif
  end subroutine eiv_diagnostics

  !-------------------------------------------------------------------------
  !> Calculate from the eigenvalue of the operator, the corresponding
  !> growth rate and frequency.
  !>
  !> \param g,f: Output of growth rate and frequency, respectively.
  !> \param lambda: Input of the eigenvalue.
  !-------------------------------------------------------------------------
  subroutine trafo_eiv_to_gf(g,f,lambda)
    use control,   only : naverage, dtim

    complex, intent(in) :: lambda
    real, intent(out)   :: g,f

    complex :: temp

    select case(mat_vec_routine)
    case(1)
      temp = log(lambda)/(naverage*dtim)
      g = real(temp)
      f = aimag(temp)
    case(2)
!      temp = log(lambda + 1.0)/(dtim)
      temp = lambda/dtim
      g = real(temp)
      f = aimag(temp)
    case default

    end select
  end subroutine trafo_eiv_to_gf

  !-------------------------------------------------------------------------
  !> Calculate from growth rate and frequency the corresponding eigenvalue
  !> of the operator.
  !>
  !> \param lambda: output parameter where to store the eigenvalue.
  !> \param g,f: Input of growth rate and frequency, respectively.
  !-------------------------------------------------------------------------
  function trafo_gf_to_eiv(g,f) result(lambda)
    use control, only : naverage, dtim
    real, intent(in) :: g,f
    complex:: lambda
    complex :: cmplx_growthrate

    select case(mat_vec_routine)
    case(1)
      cmplx_growthrate = cmplx(g, f)
      lambda = exp(naverage*dtim*cmplx_growthrate)
    case(2)
      cmplx_growthrate   = cmplx(g, f)
!      lambda = exp(dtim*cmplx_growthrate) - 1
      lambda = dtim*cmplx_growthrate
    end select
  end function trafo_gf_to_eiv

  !-------------------------------------------------------------------------
  !> Calculate from the input tolerance, the tolerance for slepc.
  !>
  !> \param tol_transformed: output parameter, the transformed tolerance.
  !> \param tolerance: input parameter, the tolerance to be transformed.
  !-------------------------------------------------------------------------
  function trafo_tolerance(tolerance) result(tol_transformed)
    use control, only : naverage, dtim
    real, intent(in) :: tolerance
    real :: tol_transformed

    select case(mat_vec_routine)
    case(1)
      tol_transformed = tolerance*dtim*naverage
    case(2)
      tol_transformed = tolerance*dtim
    end select

  end function trafo_tolerance

  !-------------------------------------------------------------------------
  !> \param ierr error return value that should be checked.
  !> \brief Checks given error value for an error and aborts if one occoured.
  !-------------------------------------------------------------------------
  subroutine check_error(ierr, where_from)
    use general, only : gkw_abort
    use global, only : int2char

    integer, intent(in) :: ierr
    character(len = *), intent(in) :: where_from

    if (0 /= ierr) then
      call gkw_abort('Error at  '//where_from&
         & //' with error code ' &
         &//int2char(ierr)//' in eigenvalue solver, aborting.')
    end if
  end subroutine check_error


#ifdef HAVE_SLEPC
  !-------------------------------------------------------------------------
  !> This function is called by Slepc to decide if the iteration should
  !> stop.
  !-------------------------------------------------------------------------
  function customEPSStoppingTest(eiv_solver,its, max_it, nconv, nev, reason)
    ! UPPERCASE USE to deliberately exclude from mkdeps script
    USE slepceps
    use control, only : max_runtime_reached
    use mpiinterface, only : mpiallreduce_or_inplace
#if defined(PETSC_USE_FORTRAN_DATATYPES)
    type(EPS) :: eiv_solver
    type(PetscInt) :: its, max_it,nconv,nev
    type(EPSConvergedReason), intent(out) :: reason
    type(PetscErrorCode) :: customEPSStoppingTest
    type(PetscErrorCode) :: ierr
#else
    EPS :: eiv_solver
    PetscInt :: its, max_it,nconv,nev
    EPSConvergedReason, intent(out) :: reason
    PetscErrorCode :: customEPSStoppingTest
    PetscErrorCode :: ierr
#endif
    logical :: stop_condition_reached

    ! do the default tests for convergence
    call EPSStoppingBasic(eiv_solver, its, max_it, nconv, nev, reason, &
       & PETSC_NULL_SCALAR, ierr)
    call check_error(ierr, 'EPSStoppingBasic')

    stop_condition_reached = max_runtime_reached(0.0) .or. stop_file_exists()
    call mpiallreduce_or_inplace(stop_condition_reached)
    if(stop_condition_reached) then
      reason = EPS_CONVERGED_USER
      call delete_stop_file()
    end if

    customEPSStoppingTest = 0

  end function customEPSStoppingTest
#endif

  function stop_file_exists()
    use control, only : stop_filename
    use mpiinterface, only : mpiwtime
    double precision, save :: t_last_check = 0.0
    logical :: stop_file_exists
    double precision :: t

    stop_file_exists = .false.

    t = mpiwtime()

    if(t - t_last_check > 10.0) then
      ! check only rarely for the file, because interaction with the
      ! filesystem is probably slow
      t_last_check = t
      inquire(file=stop_filename, exist = stop_file_exists)
      if(stop_file_exists) then
        write(*,*) 'External stop, '//stop_filename//' is present.'
        ! deleting the file here leads to a race condition.
      end if
    end if

  end function stop_file_exists

  subroutine delete_stop_file()
    use control, only : stop_filename
    use mpiinterface, only : root_processor
    if(root_processor .and. stop_file_exists()) then
      open (9, FILE = stop_filename)
      close (9, STATUS='delete')
    end if
  end subroutine delete_stop_file

end module eiv_integration
