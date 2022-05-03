!------------------------------------------------------------------------------
!>
!> DESCRIPTION
!> 
!> This diagnostic serves to output the matrix of the linear terms, as
!> computed in GKW.
!>
!------------------------------------------------------------------------------
module diagnos_matrix

  implicit none

  private
  
  public :: set_default_nml_values, bcast, check
  public :: initial_output

  !> The general on/off switch for this diagnostic
  logical, save, public :: loutput_matrix

contains

  !--------------------------------------------------------------------
  !> Set reasonable default values for the namelist items this
  !> diagnostic provides. 
  !--------------------------------------------------------------------
  subroutine set_default_nml_values()
    ! Usually the diagnostic should be switched
    ! off by default.

    loutput_matrix = .false.

  end subroutine set_default_nml_values

  !--------------------------------------------------------------------
  !> Broadcast all namelist items of this diagnostic to all processes.
  !--------------------------------------------------------------------
  subroutine bcast()
    use mpiinterface, only : mpibcast

    call mpibcast(loutput_matrix,1)
  end subroutine bcast

  !--------------------------------------------------------------------
  !> Check the diagnostic parameters and if the setup is compatible
  !> with this diagnostic.
  !--------------------------------------------------------------------
  subroutine check()

    if(.not. loutput_matrix) return

  end subroutine check

  !--------------------------------------------------------------------
  !> This routine is called at the beginning of each run (after
  !> restarts, too).
  !--------------------------------------------------------------------
  subroutine initial_output()
    use matdat, only : mat, mat_maxwll_background, mat_poisson, mat_field_diag
    use io, only : output_array, xy_fmt, binary_fmt, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use dist, only : nsolc,nf
    use mpiinterface, only : root_processor, mpireduce_sum_inplace
    real :: tmp_sum

    if(.not. loutput_matrix) return

    if(root_processor) write(*,*) 'Numbers related to the size of the problem:'
    tmp_sum = real(nsolc)
    call mpireduce_sum_inplace(tmp_sum)
    if(root_processor) then
      write(*,*) 'nsolc',int(tmp_sum)
    end if

    tmp_sum = real(nf)
    call mpireduce_sum_inplace(tmp_sum)
    if(root_processor) then
      write(*,*) 'nf_elements',int(tmp_sum)
    end if
    
    tmp_sum = real(mat%nmat)
    call mpireduce_sum_inplace(tmp_sum)
    if(root_processor) then
      write(*,*) 'nmat',int(tmp_sum)
    end if

    ! Now output only particular sections
    call output_matrix('mat', mat)
    call output_matrix('mat_maxwll_background', mat_maxwll_background)
    call output_matrix('mat_poisson', mat_poisson)
    call output_matrix('mat_field_diag', mat_field_diag)

  end subroutine initial_output

  subroutine send_matrix_to_root(matsend, matrecv, tag_start)
    use matrix_format, only : sparse_matrix, FORMAT_GKWCRS, FORMAT_GKWCOO
    use mpicomms, only : COMM_CART
    use mpiinterface, only : root_processor, send_to_root
    type(sparse_matrix), intent(inout) :: matsend, matrecv
    integer, intent(in) :: tag_start
    if(root_processor) then
      call send_to_root(matrecv%nmat,COMM_CART, tag_start+0)

      deallocate(matrecv%mat)
      allocate(matrecv%mat(matrecv%nmat))
      call send_to_root(matrecv%mat(1:matrecv%nmat), shape(matrecv%mat), COMM_CART, tag_start+1)
      
      deallocate(matrecv%jj)
      allocate(matrecv%jj(matrecv%nmat))
      call send_to_root(matrecv%jj(1:matrecv%nmat), shape(matrecv%jj), COMM_CART, tag_start+2)
      
      if(matrecv%current_format == 0) then
        deallocate(matrecv%ii)
        allocate(matrecv%ii(matrecv%nmat))
        call send_to_root(matrecv%ii(1:matrecv%nmat), shape(matrecv%ii), COMM_CART, tag_start+3)
      end if


      if(matrecv%current_format == FORMAT_GKWCRS) then
        call send_to_root(matrecv%end_row,COMM_CART, tag_start+4)
        
        deallocate(matrecv%irs)
        allocate(matrecv%irs(matrecv%end_row+1))
        call send_to_root(matrecv%irs, &
           & shape(matrecv%irs), COMM_CART, tag_start+5)
      end if

    else
      call send_to_root(matsend%nmat,COMM_CART, tag_start+0)

      call send_to_root(matsend%mat(1:matsend%nmat), &
         & shape(matsend%mat(1:matsend%nmat)), COMM_CART, tag_start+1)
      call send_to_root(matsend%jj(1:matsend%nmat), &
         & shape(matsend%jj(1:matsend%nmat)), COMM_CART, tag_start+2)
      call send_to_root(matsend%ii(1:matsend%nmat), &
         & shape(matsend%ii(1:matsend%nmat)), COMM_CART, tag_start+3)

      call send_to_root(matsend%end_row,COMM_CART, tag_start+4)
      if(matsend%current_format == FORMAT_GKWCRS) then
        call send_to_root(matsend%irs, &
           & shape(matsend%irs), COMM_CART, tag_start+5)
      end if

    end if

  end subroutine send_matrix_to_root


  subroutine output_matrix(dsetname, mat)
    use global, only : int2char_zeros
    use io, only : output_array, xy_fmt, binary_fmt, ascii_fmt, attach_metadata
    use io, only : description_key, comments_key, phys_unit_key, not_avail
    use diagnos_generic, only : attach_metadata_grid
    use general, only : gkw_warn
    use global, only : int2char
    use matrix_format, only : sparse_matrix, FORMAT_GKWCRS, FORMAT_GKWCOO
    use mpiinterface, only : root_processor, send_to_root, processor_number
    use mpiinterface, only : number_of_processors, register_tag_range, mpibarrier
    type(sparse_matrix), intent(inout) :: mat
    character(len=*), intent(in) :: dsetname
    character(len=4) :: rank_string
    integer :: rank
    
    type(sparse_matrix) :: mat_buf
    !> the range of tags to be used by this diagnostic
    integer :: tag_range_start, tag_range_end_inkl

    call register_tag_range((number_of_processors-1) * 10, tag_range_start, &
       & tag_range_end_inkl)
  
    do rank = 0, number_of_processors-1
      if(root_processor .and. rank == processor_number) then
        !write(*,*) 'root own ', mat%name, mat%nmat, mat%end_row
        mat_buf = mat
      else if (root_processor .and. rank /= processor_number) then
        !write(*,*) 'root from ',rank
        call send_matrix_to_root(mat, mat_buf, tag_range_start+rank*10)
      else if (.not.root_processor .and. rank == processor_number) then
        !write(*,*) rank, mat%name, mat%nmat, mat%end_row
        call send_matrix_to_root(mat, mat_buf, tag_range_start+rank*10)
        cycle
      end if
      rank_string = '_'//int2char_zeros(rank,3)
      
      if(root_processor) then
        call output_array(dsetname//rank_string//'_values', '/diagnostic/diagnos_matrix', &
           & mat_buf%mat(1:mat_buf%nmat), 'F', xy_fmt, binary_fmt)
        call attach_metadata(dsetname//rank_string//'_values', '/diagnostic/diagnos_matrix', &
           & phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(dsetname//rank_string//'_values', '/diagnostic/diagnos_matrix', &
           & description_key, &
           & 'the values of the nonzero elements of '//mat%name, ascii_fmt)
        call attach_metadata(dsetname//rank_string//'_values', '/diagnostic/diagnos_matrix', &
           & comments_key, &
           & 'the sparse matrix format is the compressed sparse row form&
           & (CSR or compressed row storage CRS or Yale format)', ascii_fmt)

        call output_array(dsetname//rank_string//'_jj', '/diagnostic/diagnos_matrix', &
           & real(mat_buf%jj(1:mat_buf%nmat)), 'F', xy_fmt, binary_fmt)
        call attach_metadata(dsetname//rank_string//'_jj', '/diagnostic/diagnos_matrix', &
           & phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(dsetname//rank_string//'_jj', '/diagnostic/diagnos_matrix', &
           & description_key, &
           & 'column indices', ascii_fmt)
        call attach_metadata(dsetname//rank_string//'_jj', '/diagnostic/diagnos_matrix', &
           & comments_key, &
           & not_avail, ascii_fmt)

        if(mat_buf%current_format == FORMAT_GKWCRS) then
          call output_array(dsetname//rank_string//'_irs', '/diagnostic/diagnos_matrix', &
             & real(mat_buf%irs), 'F', xy_fmt, binary_fmt)
          call attach_metadata(dsetname//rank_string//'_irs', '/diagnostic/diagnos_matrix', &
             & phys_unit_key, not_avail, ascii_fmt)
          call attach_metadata(dsetname//rank_string//'_irs', '/diagnostic/diagnos_matrix', &
             & description_key, &
             & 'the nth element of irs is the index of first element of the nth&
             & matrix row with respect to the total number of elements of the&
             & matrix. The elements of the nth matrix row are then'&
             & //dsetname//'_values(irs(n):(irs((n+1)-1)) ', ascii_fmt)
          call attach_metadata(dsetname//rank_string//'_irs', '/diagnostic/diagnos_matrix', &
             & comments_key, not_avail, ascii_fmt)
        end if

        call output_array(dsetname//rank_string//'_ii', '/diagnostic/diagnos_matrix', &
           & real(mat_buf%ii(1:mat_buf%nmat)), 'F', xy_fmt, binary_fmt)
        call attach_metadata(dsetname//rank_string//'_ii', '/diagnostic/diagnos_matrix', &
           & phys_unit_key, not_avail, ascii_fmt)
        call attach_metadata(dsetname//rank_string//'_ii', '/diagnostic/diagnos_matrix', &
           & description_key, 'row indices', ascii_fmt)
        if(mat_buf%current_format == FORMAT_GKWCOO .or. &
           & mat_buf%current_format == FORMAT_GKWCRS) then
          call attach_metadata(dsetname//rank_string//'_ii', '/diagnostic/diagnos_matrix', &
             & comments_key, not_avail, ascii_fmt)
        else
          call attach_metadata(dsetname//rank_string//'_ii', '/diagnostic/diagnos_matrix', &
             & comments_key, 'An external library provided the matrix format for this&
             & run. The generation of that datastructure may have destroyed this ii&
             & dataset.', ascii_fmt)
        end if
      end if
    end do

    call mpibarrier()

  end subroutine output_matrix
end module diagnos_matrix
