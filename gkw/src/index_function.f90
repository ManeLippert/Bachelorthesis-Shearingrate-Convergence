!-----------------------------------------------------------------------------
!> Contains the generalised index function.
!> Looks after ghost cells at processor boundaries. The data layout itself
!> is determined in dist; this module just does the indexing for a given
!> layout provided by dist.
!> The function indx(), used to obtain the index, can be called with
!> 4 (or more) arguments for fields and 6 arguments for the distribution.
!> At present, 6 is the maximum number of dimensions that can be used, but
!> it should be straightforward to extend to more if required.
!-----------------------------------------------------------------------------
module index_function

  implicit none

  private

  public :: index_init, index_set_ghostpoints, indx, dummy_index
  public :: register_offset, get_block_bounds, match_indices_with_dims, get_block_bounds_x_hack

  public :: max_dims

  integer, parameter :: max_dims = 6
  integer, parameter :: dummy_index = -9423
  integer, parameter :: max_fields = 80    !< Max number of field offsets
  integer, parameter :: max_offsets = max_dims*16

  integer, save :: index_order(max_dims)
  integer, save :: gridsize(max_dims)
  integer, save :: n_ghostpoints(max_dims,max_dims)
  integer, save :: known_map(max_offsets,max_dims+1)
  integer, save :: offset(max_offsets)
  integer, save :: field_ref(max_fields)
  integer, save :: multi_species_field(max_fields)

  integer, save :: n_offsets, n_fields, indmax
  logical, save :: register_ready = .false.

  integer, public, parameter :: IS_3D_FIELD = 2**3
  integer, public, parameter :: IS_COLL_CONS_FIELD = 2**4
  integer, public, parameter :: IS_GYROAVG_FIELD = 2**5
  integer, public, parameter :: IS_6D_DISTR = 2**6


contains


!-----------------------------------------------------------------------------
!> This routine should be called before any other routine in this module. It
!! sets the order of the dimensions in the fdisi array, calculates a
!! corresponding gridsize array to that order, sets the maximum value that
!! the index can be (based on the calculation in dist) and initialises various
!! quantities used in this module. The first 3 args for this subroutine are
!! not used but presently remain here to make the code restartable from old
!! restart files. It may be useful to allow the order to be set either via
!! input or to vary with the grid configuration. The selected index_order may
!! need to be changed after parts of the code are optimised further.
!<----------------------------------------------------------------------------

subroutine index_init(imax)

  use global,  only : id_x, id_vpar, id_mu, id_x, id_mod, id_sp, id_s
  use grid,    only : nx, nmod, nvpar, nsp, ns, nmu
  use general, only : gkw_abort
  use control, only : spectral_radius

  integer, intent(in) :: imax
  integer :: i, j

  ! The maximum index
  indmax = imax

  ! index order (to optimise mpi communications / cache efficiency)
  ! should normally reflect ordering of communicators in grid
  ! This is the default order. At present there is no way at runtime to change
  ! this. Alternatives may cause the code to run more efficiently.

  ! this ordering was previously good for linear cases, worse for nonlinear
  ! cases
!    index_order = (/ id_x, id_mod, id_s, id_vpar, id_mu, id_sp /)

  ! the order below works quite well for large collisionless nonlinear problems;
  index_order = (/ id_mod, id_x, id_s, id_vpar, id_mu, id_sp /)

  if (.not. spectral_radius) then
    ! the order below is chosen for convenience in the non-spectral case
    ! at present the gather operations in the field solve rely on this ordering
    index_order = (/ id_mod, id_s, id_x, id_vpar, id_mu, id_sp /)
  end if

  !FJC: A nice idea here is to have a runtime optimisation tests like FFTW
  ! e.g. GKW_estimate, GKW_measure, GKW_patient, GKW_exhaustive...

  ! build a correctly ordered array of grid sizes
  do i = 1, max_dims

    ! check for duplicates in the ordering
    do j = 1, max_dims
      if (i /= j .and. index_order(i) == index_order(j)) then
        call gkw_abort('index_init: duplicate index in init')
      end if
    end do

    ! Insert the correct dimension size into the current gridsize position.
    select case(index_order(i))
      case(id_vpar)
        gridsize(i) = nvpar
      case(id_s)
        gridsize(i) = ns
      case(id_sp)
        gridsize(i) = nsp
      case(id_mu)
        gridsize(i) = nmu
      case(id_x)
        gridsize(i) = nx
      case(id_mod)
        gridsize(i) = nmod
      case default
        call gkw_abort('index_init: bad index id')
    end select

  end do

  ! initialise other quantities used in this module
  n_offsets = 0
  n_fields  = 0
  n_ghostpoints(:,:) = 0
  multi_species_field(:) = 0
  field_ref(:) = 0

  ! ready to use
  register_ready = .true.

end subroutine index_init


!-----------------------------------------------------------------------------
!> Set up the number of ghost points required in each direction. The array
!! n_ghostpoints(i,j) contains the number of ghost points in the combined
!! directions i and j i.e. only for combinations of 2. If i corresponds to
!! s-direction and j corresponds to the vpar-direction, then
!! n_ghostpoints(i,i) is the number of s-points which must be communicated in
!! the s-direction and n_ghostpoints(i,j) ( == n_ghostpoints(j,i) ) is the
!! product of the number of points in s and vpar-directions that should be
!! communicated between the local processor and the processor shifted 1 in
!! both of those directions.
!<----------------------------------------------------------------------------

subroutine index_set_ghostpoints(GP_MU,GP_S,GP_X,GP_VPAR,GP_VPAR_MU,GP_VPAR_S)

  use global,  only : id_x, id_vpar, id_mu, id_x, id_mod, id_sp, id_s
  use general, only : gkw_abort

  integer, intent(in), optional :: GP_X, GP_S, GP_MU, GP_VPAR, GP_VPAR_MU, GP_VPAR_S
  integer :: i, j

  do j = 1, max_dims
    do i = 1, max_dims

      select case(index_order(i))

        case (id_mu)

          mu_outer : if (present(GP_MU)) then
            n_ghostpoints(i,i) = GP_MU
            select case (index_order(j))
              case (id_vpar)
                vpar_mu_inner : if (present(GP_VPAR_MU)) then
                  n_ghostpoints(i,j) = GP_VPAR_MU
                end if vpar_mu_inner
            end select
          end if mu_outer

        case (id_s)

          s_outer : if (present(GP_S)) then
            n_ghostpoints(i,i) = GP_S
            select case (index_order(j))
              case (id_vpar)
                vpar_s_inner : if (present(GP_VPAR_S)) then
                  n_ghostpoints(i,j) = GP_VPAR_S
                end if vpar_s_inner
            end select
          end if s_outer

        case (id_vpar)

          vpar_outer : if (present(GP_VPAR)) then
            n_ghostpoints(i,i) = GP_VPAR
            select case (index_order(j))
              case (id_s)
                s_inner : if (present(GP_VPAR_S)) then
                  n_ghostpoints(i,j) = GP_VPAR_S
                end if s_inner
              case (id_mu)
                mu_inner : if (present(GP_VPAR_MU)) then
                  n_ghostpoints(i,j) = GP_VPAR_MU
                end if mu_inner
            end select
          end if vpar_outer

        case (id_x)
        
         x_outer : if (present(GP_X)) then
            n_ghostpoints(i,i) = GP_X
            !x derivatives should not appear in combination with any other !
         end if x_outer
        
        case (id_mod,id_sp)

          ! do nothing in these cases; presently no need to index them

        case default

          call gkw_abort('index_set_ghostpoints: bad')

      end select

    end do
  end do

end subroutine index_set_ghostpoints


!-----------------------------------------------------------------------------
!> Returns an index for fdisi or a field by calling get_index(). This routine
!> re-orders the inputs beforehand as appropriate. The first argument passed
!> to this routine is an identifier, which checks if the correct number of
!> subsequent arguements has been provided, and determines their meaning in
!> terms of grid indices.
!>
!> The returned index refers to the single solution vector storing everything
!> 
!> Example calls:
!>  fdis = fdisi(indx(ifdis,imod,ix,i,j,k,is))
!>  phi  = fdisi(indx(iphi,imod,ix,i))
!>
!> where standard conventions used for the indices are
!>    imod, ix = mode position
!>           i = s position
!>           j = mu position
!>           k = vpar position
!>          is = species number
!>
!> The dimension indices always appear in the order above as for fdisi
!> but for lower dimensional fields those not used are ommitted.
!> 
!> The identifiers (defined in dist) should have certain bits set:
!>   3D fields: iphi, iapar, ibpar   (imod,ix,i)        : bit IS_3D_FIELD
!>   4D momentum change field: i_mom (imod,ix,i,is)     : bit IS_COLL_CONS_FIELD
!>   5D gyro-averages: iphi_ga, etc  (imod,ix,i,j,is)   : bit IS_GYROAVG_FIELD
!>   6D distributions: ifdis         (imod,ix,i,j,k,is) : bit IS_6D_DISTR
!>
!APS: These can be used to check the registration of each component in dist,
!APS: before the index function is called.
!-----------------------------------------------------------------------------

function indx(i_field,I1,I2,I3,I4,I5,I6)

  use global,  only : id_sp, id_s, id_vpar, id_mu, id_mod, id_x
  use general, only : gkw_abort

  integer,           intent(in) :: i_field
  integer, optional, intent(in) :: I1, I2, I3, I4, I5, I6
  integer                       :: indx

  integer :: ind_array(max_dims)
  integer :: i

  ! insert dummies to start with (so that redundant indices can be identified)
  ind_array(:) = dummy_index

  ! Check if the first argument is a valid identifier, and if so, also check
  ! that the correct number of arguements have been passed.

  if (iand(i_field,IS_6D_DISTR) /= 0) then

    !
    ! full distribution function
    !   location = indx(identifier,imod,ix,i,j,k,is)
    !

    if (present(I1) .and. present(I2) .and. present(I3) .and. present(I4)    &
        &           .and. present(I5) .and. present(I6)) then
      ! called with the correct number of inputs; order them
      do i = 1, max_dims
        ! The index_order array determines the order in which the indices are
        ! stored in the array of indices for indx evaluation. This is *not*
        ! related to the order in which the arguments appear in the call to
        ! the function indx(), which is reflected in the fixed mapping below
        ! between ids and argument number.
        select case(index_order(i))
        case (id_mod)
          ind_array(i) = I1
        case (id_x)
          ind_array(i) = I2
        case (id_s)
          ind_array(i) = I3
        case (id_mu)
          ind_array(i) = I4
        case (id_vpar)
          ind_array(i) = I5
        case (id_sp)
          ind_array(i) = I6
        case default
          call gkw_abort('bad index_6d') ! will never happen
        end select
      end do
    else
      write(*,*) 'indx: expected 6 indicies for identifier', i_field
      call gkw_abort('indx: called with wrong number of indicies')
    end if

  else if (iand(i_field,IS_GYROAVG_FIELD)/=0) then

    !
    ! gyro average
    !   location = indx(identifier,imod,ix,i,j,is)
    !

    if (present(I1) .and. present(I2) .and. present(I3) .and. present(I4)    &
        &           .and. present(I5) .and. (.not. present(I6))) then
      ! called with the correct number of inputs; order them
      do i = 1, max_dims
        select case(index_order(i))
        case (id_mod)
          ind_array(i) = I1
        case (id_x)
          ind_array(i) = I2
        case (id_s)
          ind_array(i) = I3
        case (id_mu)
          ind_array(i) = I4
        case (id_sp)
          ind_array(i) = I5
        case (id_vpar)
          ! do nothing
        case default
          call gkw_abort('bad index_5d') ! will never happen
        end select
      end do
    else
      write(*,*) 'indx: expected 5 indicies for identifier', i_field
      call gkw_abort('indx: called with wrong number of indicies')
    end if

  else if (iand(i_field,IS_COLL_CONS_FIELD) /= 0) then

    !
    ! momentum conservation
    !   location = indx(identifier,imod,ix,i,is)
    !

    if (present(I1) .and. present(I2) .and. present(I3) .and. present(I4)    &
        &           .and. (.not. present(I5)) .and. (.not. present(I6))) then
      ! called with the correct number of inputs; order them
      do i = 1, max_dims
        select case(index_order(i))
        case (id_mod)
          ind_array(i) = I1
        case (id_x)
          ind_array(i) = I2
        case (id_s)
          ind_array(i) = I3
        case (id_sp)
          ind_array(i) = I4
        case (id_mu,id_vpar)
          ! do nothing
        case default
          call gkw_abort('bad index_4d') ! will never happen
        end select
      end do
    else
      write(*,*) 'indx: expected 4 indicies for identifier', i_field
      call gkw_abort('indx: called with wrong number of indicies')
    end if

  else if (iand(i_field,IS_3D_FIELD)/=0) then

    !
    ! regular field
    !   location = indx(identifier,imod,ix,i)
    !

    if (         present(I1)  .and.        present(I2)  .and.                &
        &        present(I3)  .and. (.not. present(I4)) .and.                &
        & (.not. present(I5)) .and. (.not. present(I6))) then
      ! called with the correct number of inputs; order them
      do i = 1, max_dims
        select case(index_order(i))
        case (id_mod)
          ind_array(i) = I1
        case (id_x)
          ind_array(i) = I2
        case (id_s)
          ind_array(i) = I3
        case (id_mu,id_vpar,id_sp)
          ! do nothing
        case default
          call gkw_abort('bad index_3d') ! will never happen
        end select
      end do
    else
      write(*,*) 'indx: expected 3 indicies for identifier', i_field
      call gkw_abort('indx: called with wrong number of indicies')
    end if

  else
  
    ! unknown; abort
    write(*,*) 'indx: bad identifier ', i_field
    call gkw_abort('indx: called with unknown identifier')

  end if
    
  ! abort if the switch is unknown - should have been registered
  if (.not. any(field_ref(1:n_fields) == i_field)) then
    write (*,*) 'indx: unknown switch', i_field
    write (*,*) 'known field refs:', field_ref(1:n_fields)
    call gkw_abort('index_other: unregistered switch')
  end if

  ! call the general index function
  indx = get_index(ind_array,i_field)

end function indx


!-----------------------------------------------------------------------------
!> Routine that works out the index. ind_array contains 6 integers
!! corresponding to the grid points in the various directions. When all the
!! integers are in the normal range, from 1 to the number of local grid points
!! in the direction, the index is
!!  ind_array(1) + (ind_array(2) - 1)*gridsize(1) + ...
!!               + (ind_array(6) - 1)*gridsize(1)*gridsize(2)*...*gridsize(5),
!! where gridsize is the corresponding array of local grid sizes. This is the
!! case for the full local part of the solution, fdisi(1:nf).
!!
!! If the ind_array values are outside the normal range by a few, the call may
!! refer to a point of the solution on another processor. If the input refers
!! to the solution on the next processor in some direction, then it will be
!! larger than the maximum gridsize value in that direction. When the size of
!! the grid in that direction is subtracted from that value, it should then be
!! between 1 and the maximum number of `ghost' points that we would wish to
!! reference. If the input value is too small, it can be shifted up into this
!! same range. To generate an index for the bit of the array that contains
!! the points from a particular adjacent processor, the gridsize entry for
!! that direction can be reduced to the number of ghost points in that
!! direction, then the above can be used with the shifted ind_array value.
!! The block containing the solution from another processor is offset from the
!! start of fdisi by some amount which is registered from dist. The pattern
!! of shifts in ind_array is stored in the remap array, then used to look up
!! the shift provided by dist.
!!
!! In the case of fields, which only span 3 or so of the dimensions, a dummy
!! index is located in unused dimensions when this routine is called. The
!! shifts and offsets can be obtained in the same way as the distribution
!! function. The dummy dimensions are ignored for the purpose of calculating
!! an index value. The switch is stored in an extra slot in the remap array.
!<----------------------------------------------------------------------------

function get_index(index_arr,i_field)

  use general, only : gkw_abort

  integer, intent(in) :: index_arr(max_dims)
  integer, intent(in) :: i_field
  integer :: get_index

  integer :: fact(max_dims), length(max_dims), ind_array(max_dims)
  integer :: map(max_dims+1)
  integer :: i, constant, blck, ioffset, j, outside_range_count
  integer :: jloc

  ! initialise the length array to the array of gridsizes
  length = gridsize

  ind_array = index_arr

  ! initialise the map and factor arrays
  map(:) = 0 ; fact(:) = 0

  !
  ! perform a check and shift (if necessary) on every input dimension
  !

  ! First count how many indices are out of range. The maximum is 2.
  ! because we have mixed derivatives of no more than second order,
  ! involving 2 different directions.
  outside_range_count = 0
  do i = 1, max_dims
    if ((ind_array(i) < 1 .or. ind_array(i) > length(i)) .and.               &
        &  ind_array(i) /= dummy_index) then
      outside_range_count = outside_range_count + 1
    end if
  end do
  if (outside_range_count > 2) then
    ! then have coded wrong
    call gkw_abort('get_index: there are presently valid calls with more'//  &
                  &'than 2 indices out of range')
  end if

  check_all_dims : do i = 1, max_dims

    ! Ignore dummy entries; this array may contain dummies if referencing
    ! a field where some entries are meaningless.
    not_dummy : if (ind_array(i) /= dummy_index) then

      ! set a default zero for the map
      map(i) = 0

      ! Check if the input value is less than 1
      too_small : if (ind_array(i) < 1) then
        ! If there are ghost points in this direction, it could refer to a
        ! point on the previous processor. Shift the value and check.
        processor_below : if (n_ghostpoints(i,i) > 0) then
          if (outside_range_count == 1) then
            ind_array(i) = ind_array(i) + n_ghostpoints(i,i)
            length(i)    = n_ghostpoints(i,i)
            ! Check that the shifted value is now in the acceptable range for
            ! a ghost point i.e. between 1 and n_ghostpoints(i,i).
            if (ind_array(i) < 1 .or. ind_array(i) > n_ghostpoints(i,i)) then
              call gkw_abort('indx main: bad input')
            end if
          else if (outside_range_count == 2) then
            ! find the other index which is shifted
            jloc = 0
            do j = 1, max_dims
              if (i /= j) then
                if (n_ghostpoints(i,j) /= 0) then
                  jloc = j
                end if
              end if
            end do
            if (jloc == 0) call gkw_abort('get_index: bad jloc below')
            ind_array(i) = ind_array(i) + n_ghostpoints(i,jloc)
            length(i)    = n_ghostpoints(i,jloc)
            ! Check that the shifted value is now in the acceptable range for
            ! a ghost point i.e. between 1 and n_ghostpoints(i,i).
            if (ind_array(i) < 1 .or. ind_array(i) > n_ghostpoints(i,jloc))  &
                & call gkw_abort('indx main: bad input')

          end if
          ! Make a note that this point is on the previous processor.
          map(i) = -1
        else ! a bad input value
          call gkw_abort('indx_main: i < 1')
        end if processor_below

      end if too_small

      ! Check if the input value is greater than the number of points in the
      ! local grid for the present direction.
      too_big : if (ind_array(i) > length(i)) then
        ! If there are ghost points in this direction, the point may be on the
        ! next processor. Shift the value into the right range.
        processor_above : if (n_ghostpoints(i,i) > 0) then
          if (outside_range_count == 1) then
            ind_array(i) = ind_array(i) - length(i)
            length(i)    = n_ghostpoints(i,i)
            ! Check that the shifted value is now in the acceptable range for
            ! a ghost point i.e. between 1 and n_ghostpoints(i,i).
            if (ind_array(i) < 1 .or. ind_array(i) > n_ghostpoints(i,i)) then
              call gkw_abort('indx main: bad input')
            end if
          else if (outside_range_count == 2) then
            ! find the other index which is shifted
            jloc = 0
            do j = 1, max_dims
              if (i /= j) then
                if (n_ghostpoints(i,j) /= 0) then
                  jloc = j
                end if
              end if
            end do
            if (jloc == 0) call gkw_abort('get_index: bad jloc below')
            ind_array(i) = ind_array(i) - length(i)
            length(i)    = n_ghostpoints(i,jloc)
            ! Check that the shifted value is now in the acceptable range for
            ! a ghost point i.e. between 1 and n_ghostpoints(i,i).
            if (ind_array(i) < 1 .or. ind_array(i) > n_ghostpoints(i,jloc))  &
                & call gkw_abort('indx main: bad input')

          end if
          ! Make a note that this point is on the next processor
          map(i) = 1
        else ! a bad input value
          write(*,*) ind_array
          write(*,*) length 
          call gkw_abort('indx_main: i > nmax')
        end if processor_above

      end if too_big

    else

      ! put a dummy in the length array
      length(i) = dummy_index

    end if not_dummy

  end do check_all_dims

  !
  ! compact the arrays e.g. (12,4,dummy,3,dummy,dummy) -> (12,4,3,0,0,0)
  ! This allows the different types of blocks to be dealt with in the same
  ! way.
  !

  call compact_array(length)
  call compact_array(ind_array)

  !
  ! work out the array (fact) and integer (constant) to calculate indx
  !
  ! for the above Example length=(12,4,3,0,0,0) this will calculate
  ! constant = 1 -1 -1*12 - 1*12*4 - 1*12*4*3
  ! fact = (1, 1*12, 1*12*4, 1*12*4*3)

!APS: It should be possible to pre-calculate these constants for each embedded
!APS: array in dist; they are (almost) independent of the input indices. The
!APS: simplest thing to do would be to store the constants for a known "length"
!APS: array and re-use them on the next call. It would be better to do it via
!APS: register_offset() if possible (which it should be, since this routine
!APS: generates the "map" array that is used to generate the pattern used
!APS: to recall the offset).

  ! initial values
  blck = 1
  constant = 1

  set_blocksizes : do i = 1, max_dims

    ! Exit if the length is zero, otherwise the bit to add on will be wrong.
    if (length(i) == 0) exit set_blocksizes

    ! set the prefactors and update the additional bit
    fact(i)  = blck
    constant = constant - blck

    ! update the block size
    blck = blck*length(i)

  end do set_blocksizes

  !
  ! Calculate the index, using the map array to obtain the pre-registered
  ! block offset.
  !

  ! store any i_field for a `field' in an additional slot in map
  map(max_dims+1) = i_field

  ! get the offset
  ioffset = index_offset(map)

  get_index = dot_product(fact,ind_array) + constant + ioffset

  ! check for out of bounds
  if (get_index < 1 .or. get_index > indmax) then
    write (*,*) 'get_index: (indx, i_field) =', get_index, i_field
    call gkw_abort('get_index: bad indx value')
  end if

end function get_index


!-----------------------------------------------------------------------------
!> Obtain the offset in fdisi by making a comparison of the shifts and
!> switches recorded in the map array with those in a table of patterns and
!> offsets.
!-----------------------------------------------------------------------------

function index_offset(map)
  use general, only : gkw_abort
  integer, intent(in) :: map(max_dims+1)
  integer :: index_offset
  integer :: i

  ! loop over the known offsets
  known_patterns : do i = 1, n_offsets
    ! check if the generated pattern matches the list entry
    if (all(map == known_map(i,:))) then
      ! return the offset if the pattern matches
      index_offset = offset(i)
      return
    end if
  end do known_patterns

  ! The pattern is not in the list, so no offset can be returned.
  write (*,*) 'map:       ',map(:)
  call gkw_abort('index_offset: unregistered offset')

end function index_offset


!-----------------------------------------------------------------------------
!> Record the offset together with call pattern. From dist, a call to this
!> routine is required from each block. The offset/pattern combination is
!> then used to get the right offset when the index function is called. The
!> default pattern for any direction is 0 (a zero) if it is not present (all
!> directions are optional). The combination is stored in known_map.
!-----------------------------------------------------------------------------

subroutine register_offset(IMOD,IX,IS,IMU,IVPAR,ISP,ioffset,IFIELD,NSP)

  use global,  only : id_x, id_s, id_sp, id_vpar, id_mu, id_mod
  use general, only : gkw_abort

  integer, optional, intent(in) :: IS, ISP, IMU, IVPAR, IX, IMOD, NSP, IFIELD
  integer, intent(in) :: ioffset

  integer :: map(max_dims+1)
  integer :: i, i_vpar, i_mu, i_s, i_sp, i_x, i_mod

  if (.not. register_ready) call gkw_abort('register_offset: not ready')

  map(:) = 0
  ! all zero defaults
  i_mod = 0 ; i_s = 0 ; i_sp = 0 ; i_vpar = 0 ; i_x = 0 ; i_mu = 0

  ! change patterns for items present
  if (present(IMOD))   i_mod = IMOD
  if (present(IS))       i_s = IS
  if (present(ISP))     i_sp = ISP
  if (present(IMU))     i_mu = IMU
  if (present(IVPAR)) i_vpar = IVPAR
  if (present(IX))       i_x = IX

  ! for the fields, make a fields reference list and store the switch.
  if (present(IFIELD)) then
    n_fields = n_fields + 1
    field_ref(n_fields) = IFIELD
    if (present(NSP))then
      multi_species_field(n_fields) = NSP
    end if
    ! the last map value contains the switch
    map(max_dims+1) = IFIELD
  end if

  ! create the pattern
  do i = 1, max_dims
    select case(index_order(i))
      case (id_sp)
        map(i) = i_sp
      case (id_s)
        map(i) = i_s
      case (id_x)
        map(i) = i_x
      case (id_mod)
        map(i) = i_mod
      case (id_mu)
        map(i) = i_mu
      case (id_vpar)
        map(i) = i_vpar
      case default
        write (*,*) index_order
        write (*,*) map
        call gkw_abort('bad register offset') ! will never happen
    end select
  end do

  ! update the list of offsets/patterns
  n_offsets = n_offsets + 1
  known_map(n_offsets,:) = map(:)
  offset(n_offsets) = ioffset

end subroutine register_offset


!-----------------------------------------------------------------------------
!> Take out the dummy entries and shift all real entries to the left; put a
!> zero in the remaining postitions. Note that the grab_dummies loop must be
!> greedy and keep removing dummies till they are gone.
!-----------------------------------------------------------------------------

subroutine compact_array(iarray)

  integer, intent(inout) :: iarray(max_dims)
  integer :: i, j

  do i = 1, max_dims
    grab_dummies : do
      if (iarray(i) == dummy_index) then
        do j = i, max_dims-1
          iarray(j) = iarray(j+1)
        end do
        iarray(max_dims) = 0
      else
        exit grab_dummies
      end if
    end do grab_dummies
  end do

end subroutine compact_array


!-----------------------------------------------------------------------------
!>
!> If called without any optional arguments, this routine returns something like
!>  starts(:) = 1
!>  ends(:) = (/ nmod, nx, ns, nvpar, nmu, nsp/)
!> but ordered according to index_order.
!>
!> If called with is_distr_or_field argument, the return values are
!> either
!>  ends(:) = (/ 1   , nx, ns, 1, 1, 1/)
!> or (for IS_GYROAVG_FIELD)
!>  ends(:) = (/ nmod, nx, ns, 1, nmu, nsp/)
!>
!> If called with any _NEXT argument, the corresponding entry in starts(:)
!> becomes, for example
!>   = nvpar - GPX_NEXT + 1
!>
!> If called with any _PREV argument, the corresponding entry in ends(:)
!> becomes, for example
!>   = GPX_PREV
!> For
!> 
!>
!>
!>
!-----------------------------------------------------------------------------

subroutine get_block_bounds(starts,ends,is_distr_or_field,&
                       & GPS_NEXT,GPS_PREV,GPVPAR_NEXT,GPVPAR_PREV,          &
                       & GPMU_PREV,GPMU_NEXT,GPX_NEXT, GPX_PREV)

  use global,  only : id_sp, id_s, id_vpar, id_mu, id_mod, id_x
  use grid,    only : nx, nmod, nvpar, ns, nmu, nsp
  use general, only : gkw_abort

  integer, intent(in), optional :: GPS_NEXT, GPS_PREV, GPVPAR_PREV, GPMU_NEXT
  integer, intent(in), optional :: GPVPAR_NEXT, GPMU_PREV, GPX_PREV, GPX_NEXT
  !> this should be either IS_6D_DISTR, IS_3D_FIELD or IS_GYROAVG_FIELD
  !> N.B. we have not considered communicated mom conserve types
  integer, intent(in) :: is_distr_or_field
  integer, intent(out) :: starts(max_dims), ends(max_dims)

  integer :: i

  starts(:) = 1
  ends(:)   = 1

  do i = 1, max_dims
    ! build the default length array
    select case(index_order(i))
      case(id_vpar)
        if (present(GPVPAR_PREV)) then
          ends(i) = GPVPAR_PREV
        else
          if (is_distr_or_field == IS_6D_DISTR) then
            ends(i) = nvpar
          else
            ends(i) = 1
            ! this should never be asked for
            !stop 'wrong'
          end if
        end if
        if (present(GPVPAR_NEXT)) then
          starts(i) = nvpar - GPVPAR_NEXT + 1
        end if
      case(id_s)
        if (present(GPS_PREV)) then
          ends(i) = GPS_PREV
        else
          ends(i) = ns
        end if
        if (present(GPS_NEXT)) then
          starts(i) = ns - GPS_NEXT + 1
        end if
      case(id_sp)
        select case(is_distr_or_field) 
        case(IS_3D_FIELD)
          ends(i) = 1
        case(IS_GYROAVG_FIELD)
          ends(i) = nsp
        case(IS_6D_DISTR)
          ends(i) = nsp
        end select
      case(id_mu)
        if (present(GPMU_PREV)) then
          ends(i) = GPMU_PREV
        else
          select case(is_distr_or_field)
          case(IS_3D_FIELD)
            ends(i) = 1
          case(IS_GYROAVG_FIELD)
            ends(i) = nmu
          case(IS_6D_DISTR)
            ends(i) = nmu
          end select
        end if
        if (present(GPMU_NEXT)) then
          starts(i) = nmu - GPMU_NEXT + 1
        end if
      case(id_x)
        if (present(GPX_PREV)) then
          ends(i) = GPX_PREV
        else
          ends(i) = nx
        end if
        if (present(GPX_NEXT)) then
          starts(i) = nx - GPX_NEXT + 1
        end if
      case(id_mod)
        ends(i) = nmod
      case default
        call gkw_abort('get_block_bounds_: bad index id')
    end select
  end do

end subroutine get_block_bounds


!-----------------------------------------------------------------------------
!> A hack to create a ghost recv datatype for non contiguous blocks.
!> This routine overwrites relevant starts and ends after calling get_block_bounds with
!> no ghost point argument.  This functionality could be merged with that 
!> routine if it is made a little more general.
!-----------------------------------------------------------------------------
subroutine get_block_bounds_x_hack(starts,ends,gp_pm)

  use global,  only : id_x
  use grid,    only : nx 

  integer, intent(in) ::  gp_pm 
  integer, intent(inout) :: starts(max_dims), ends(max_dims)

  integer :: i

   do i = 1, max_dims
     if (index_order(i) == id_x) then
       if (gp_pm > 0) then
         starts(i) = nx + 1
         ends(i)   = nx + gp_pm
       else if (gp_pm < 0) then
         starts(i) = gp_pm + 1
         ends (i)  = 0
       end if
     end if     
   end do

end subroutine get_block_bounds_x_hack

!-----------------------------------------------------------------------------
!> Re-orders the re-verted index before calling
!-----------------------------------------------------------------------------

subroutine match_indices_with_dims(i1,i2,i3,i4,i5,i6,i_mod,i_x,i_s,i_mu,i_vpar,i_sp)

  use global,  only : id_sp, id_s, id_vpar, id_mu, id_mod, id_x
  use general, only : gkw_abort

  integer, intent(in)  :: i1, i2, i3, i4, i5, i6
  integer, intent(out) :: i_mod, i_x, i_s, i_mu, i_vpar, i_sp

  integer :: ind_in(max_dims)
  integer :: i

  ind_in(1:6) = (/ i1, i2, i3, i4, i5, i6 /)

  do i = 1, max_dims
    select case(index_order(i))
      case(id_vpar)
        i_vpar = ind_in(i)
      case(id_s)
        i_s    = ind_in(i)
      case(id_mu)
        i_mu   = ind_in(i)
      case(id_x)
        i_x    = ind_in(i)
      case(id_sp)
        i_sp   = ind_in(i)
      case(id_mod)
        i_mod  = ind_in(i)
      case default
        ! should not happen
        call gkw_abort('match_indices_with_dims: unknown id')
    end select
  end do

end subroutine match_indices_with_dims


end module index_function
