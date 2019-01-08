subroutine dump_mat(factidx, n, ne, row, col, a)
   implicit none

   integer, intent(in) :: factidx, n, ne
   integer, intent(in) :: row(ne), col(ne)
   double precision, intent(in) :: a(ne)

   character(len=80) :: filename
   character(len=72) :: title
   character(len=8) :: key
   integer :: icntl(10), info(5)

   integer :: i
   integer, allocatable :: ip(:), ind(:), iw(:)
   double precision, allocatable :: anew(:)

   icntl(1) = 12 ! unit number
   icntl(2) = 17 ! output precision
   icntl(3) = 1 ! values supplied (not just pattern)
   icntl(4) = 0 ! symmetric matrix
   icntl(5) = 1 ! coordinate form

   write(key, "(i8)") factidx
   key = adjustl(key)
   title = "IPOPT MA57 dump factorization number " // trim(key)
   filename = "matrix." // trim(key) // ".rb"
   open(icntl(1), file=filename, status="replace")

   allocate(ip(2*n+1), ind(2*ne), iw(n+1), anew(ne))
   ind(1:ne) = row(:)
   ind(ne+1:2*ne) = col(:)
   anew(:) = a(:)

   call mc54ad(icntl, title, key, n, n, ne, ip, ind, anew, iw, info)

   close(icntl(1))
end subroutine dump_mat

subroutine dump_mat_csc(factidx, n, ptr, row, a)
   implicit none

   integer, intent(in) :: factidx, n
   integer, intent(in) :: ptr(n+1), row(ptr(n+1)-1)
   double precision, intent(in) :: a(ptr(n+1)-1)

   character(len=80) :: filename
   character(len=72) :: title
   character(len=8) :: key
   integer :: icntl(10), info(5)

   integer :: i, ne
   integer, allocatable :: ip(:), ind(:), iw(:)
   double precision, allocatable :: anew(:)

   ne = ptr(n+1)-1

   icntl(1) = 12 ! unit number
   icntl(2) = 17 ! output precision
   icntl(3) = 1 ! values supplied (not just pattern)
   icntl(4) = 0 ! symmetric matrix
   icntl(5) = 0 ! CSC form

   write(key, "(i8)") factidx
   key = adjustl(key)
   title = "IPOPT MA86 dump factorization number " // trim(key)
   filename = "matrix." // trim(key) // ".rb"
   open(icntl(1), file=filename, status="replace")

   allocate(ip(2*n+1), ind(2*ne), iw(n+1), anew(ne))
   ip(1:n+1) = ptr(1:n+1)
   ind(1:ne) = row(:)
   anew(:) = a(:)

   call mc54ad(icntl, title, key, n, n, ne, ip, ind, anew, iw, info)

   close(icntl(1))
end subroutine dump_mat_csc
