!-----------------------------------------------------------------------
!> module that contains routines for converting matrix formats
!-----------------------------------------------------------------------
module matconv

  use global, only : iumf, rumf
  use general, only : gkw_abort

implicit none

#if defined(umfpack)

public :: crstoreal, crstoreal_umf

private

contains
!-----------------------------------------------------------------
!>     this routine converts the complex row storage format to
!>     a real row storage format
!>
!>     np     number of columns ("rows") in the complex matrix
!>     nm     maximum dimension of mat and jam
!>     iam    integer array that points to the beginning of 
!>            the next row in the mat array
!>     jam    index array that points to the column entry
!>     mat    the matrix elements
!>     nmr    maximum dimension of a and ja
!>  
!>     output 
!>     ia     same as iam but for the real matrix (twice the size)
!>     ja     same as jam but for the real matrix (4 times in size)
!>     a      same as mat but for the real matrix (4 times in size)
!-------------------------------------------------------------------
!> This routine could be further optimised to deal with real 
!> entries in the matrix (which give zeros on some of the 
!> diagonals
subroutine crstoreal(np,nm,iam,jam,mat,nmr,ia,ja,a)
  use global, only : gkw_a_equal_b_accuracy
integer          :: np, nm, iam(np+1), jam(nm) 
integer(kind=iumf)        :: nmr, itel 
integer(kind=iumf)        :: ia(2*np+1), ja(nmr)
complex          :: mat(nm) 
real(kind=rumf)  :: a(nmr)

integer(kind=iumf) i, j, istart

! start the conversion 
do i = 1, np 
  ia(2*i) = 2*(iam(i+1)-iam(i))
  ia(2*i+1) = 2*(iam(i+1)-iam(i))
end do 
ia(1) = 1
do i = 2, 2*np+1
  ia(i) = ia(i) + ia(i-1)
end do 
itel = 1
do i = 1, np
  do j = iam(i), iam(i+1)-1
    ja(itel) = 2*jam(j)-1
    a(itel) = real(mat(j))
    itel = itel+1
    ja(itel) = 2*jam(j)
    a(itel) = -aimag(mat(j))
    itel = itel+1
  end do 
  do j = iam(i), iam(i+1)-1
    ja(itel) = 2*jam(j)-1
    a(itel) = aimag(mat(j))
    itel = itel+1
    ja(itel) = 2*jam(j)
    a(itel) = real(mat(j))
    itel = itel+1
  end do 
end do 

if (.true.) then 
! compress the newly generated matrix (many zeros might 
! have appeared if the matrix elements were real) 
itel = 0
ia(1) = 1 
istart = 1
do i = 1, 2*np 
  do j = istart, ia(i+1)-1
    if (.not. gkw_a_equal_b_accuracy(real(a(j)), 0.0)) then
      itel = itel + 1 
      a(itel) = a(j) 
      ja(itel) = ja(j) 
    endif 
  end do 
  istart = ia(i+1)
  ia(i+1) = itel + 1 
end do 

endif 

end subroutine crstoreal 

!-----------------------------------------------------------------
!>     this routine converts the complex row storage format to
!>     a real row storage format and also makes any type 
!>     conversions necessary for umfpack
!>
!>     input (in the precision of GKW)
!>     np     number of colums in the complex matrix
!>     nm     maximum dimension of mat and jam
!>    iam    integer array that points to the beginning of 
!>            the next row in the mat array
!>     ii     index array that points to the column entery
!>     mat    the matrix elements
!>     nmr    maximum dimension of a and ja
!>  
!>     outputs (always in double precision, uses umf kinds) 
!>     apk    same as iam but for the real matrix (twice the size)
!>    aik    same as jam but for the real matrix (4 times in size)
!>   axk    same as mat but for the real matrix (4 times in size)
!------------------------------------------------------------------
!> This routine could be further optimised to deal with real 
!> entries in the matrix (which give zeros on some of the 
!> diagonals
subroutine crstoreal_umf(np,iam,mat,nmr,apk,aik,axk)
  use global, only : gkw_a_equal_b_accuracy
  use matrix_format, only : sparse_matrix
  integer,intent(in) :: np, iam(np+1)
  integer(kind=iumf),intent(in) :: nmr
  integer(kind=iumf),intent(out) :: apk(2*np+1), aik(nmr)
  type(sparse_matrix),intent(in) :: mat
  real(kind=rumf),intent(out) :: axk(nmr)

  integer(kind=iumf) :: itel
  integer(kind=rumf) :: i, j, istart

if (nmr < 0) call gkw_abort('error calling crstoreal_umf')
if (2*np < 0) call gkw_abort('error calling crstoreal_umf')

! start the conversion 
do i = 1, np 
  apk(2*i) = 2*(iam(i+1)-iam(i))
  apk(2*i+1) = 2*(iam(i+1)-iam(i))
end do 
apk(1) = 1
do i = 2, 2*np+1
  apk(i) = apk(i) + apk(i-1)
end do 

itel = 1
do j = 1, np
  do i = iam(j), iam(j+1)-1
    aik(itel) = 2*mat%ii(i)-1
    axk(itel) = real(mat%mat(i),rumf)
    itel = itel+1
    aik(itel) = 2*mat%ii(i)
    axk(itel) = real(aimag(mat%mat(i)),rumf)
    itel = itel+1
    if (itel+2 > huge(itel)-10) call gkw_abort('integer overflow')
  end do 
  do i = iam(j), iam(j+1)-1
    aik(itel) = 2*mat%ii(i)-1
    axk(itel) = -real(aimag(mat%mat(i)),rumf)
    itel = itel+1
    aik(itel) = 2*mat%ii(i)
    axk(itel) = real(mat%mat(i),rumf)
    itel = itel+1
    if (itel > huge(itel)-10) call gkw_abort('integer overflow')
  end do 
end do 

if (.true.) then 
! compress the newly generated matrix (many zeros might 
! have appeared if the matrix elements were real) 
itel = 0
apk(1) = 1 
istart = 1
do j = 1, 2*np 
  do i = istart, apk(j+1)-1
    if (.not. gkw_a_equal_b_accuracy(real(axk(i)), 0.0)) then
      itel = itel + 1 
      axk(itel) = axk(i) 
      aik(itel) = aik(i) 
    endif 
  end do 
  istart = apk(j+1)
  apk(j+1) = itel + 1 
end do 

endif 

end subroutine crstoreal_umf

#endif

end module matconv
