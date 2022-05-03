!-------------------------------------------------------------------------------
!> Module that defines and stores structures used in GKW.
!-------------------------------------------------------------------------------
module structures

  implicit none

  public

  type :: matrix_element
    !> IMPORTANT NOTE: The MPI datatype in module mpiinterface relies on the
    !> precise order of these struct members. Do NOT CHANGE it without
    !> correcting the construction of that MPI datatype.


    !> coordinates specifying the "row index" of the matrix element
    integer :: imod   !< the toroidal mode number 
    integer :: ix     !< the radial grid point 
    integer :: i      !< the grid point along the field 
    integer :: j      !< the grid point in mu 
    integer :: k      !< the grid point in vpar 
    integer :: is     !< the grid point in the species

    !> coordinates specifying the "column index" of the matrix element
    integer :: imloc  !< the desired toroidal mode number 
    integer :: ixloc  !< the desired location along the radial direction
    integer :: iloc   !< the desired location along the s-direction
    integer :: jloc   !< the desired location along the mu-direction
    integer :: kloc   !< the desired location along the vpar-direction
    integer :: isloc  !< the desired species index 
    ! The following components have default values which must 
    ! be overwritten in the prototype element before register_term 
    ! is called; matdat checks against these defaults.
    ! These numbers must NOT overlap with identifier ranges in dist.
    integer :: itype = -300
                      !< type of the iih index (ifdis,iphi,iapar,ibpar,i_mom)
    integer :: itloc = -200
                      !< type of the jjh index (ifdis,iphi,iapar,ibpar,i_mom)
    integer :: ideriv = -100
                      !< order of derivative (0,1,2,4) for timestep estimate
                      !< (NOT the order of the differential scheme)

    ! labels for special elements
    logical :: outflow = .false.  !< outflow elements for energetics diagnostic

    complex :: val    !< the value of the element
    !> the name of the term that put the element
    character(len=64) :: term =  ''
  end type matrix_element

  !> structure for storing MPI persistent communication information.
  !> this is used in exp integration for communication of fdisi ghost cells.
  !> allocatable attribute is not allowed in fortan95 inside structures.
  type :: ghostcomm
    integer :: nreq = 0                            !< number of requests total
    integer :: irequest = 0                        !< number of requests counter
    integer, pointer, dimension (:) :: preq        !< status returned by mpiwait
    !integer, allocatable, dimension (:,:) :: pstat !< requests labels
    logical :: running = .false.                   !< true only when running
    logical :: initialised = .false.               !< true only when initialised
  end type ghostcomm

end module structures
