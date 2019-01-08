! F90 dependencies
! COPYRIGHT (c) 2006 Council for the Central Laboratory
!               of the Research Councils
! Original date 21 August 2006. Version 1.0.0.
! Version 1.1.0. 11 October 2007. Optional argument partial added to three
!                subroutines.
! 2 November 2007. Version 1.1.0.
!                  Optional argument partial added to three subroutines.
!                  Coding of ma54_to_block and ma54_from_block improved.
! 29 April 2008. Version 1.2.0.
!                Uses long integers internally so that n can be larger.
! 13 August 2008. Version 1.3.0.
!               Optional argument partial ignored.
!               Trailing part rearranged only temporarily while it is updated.
!               OpenMp support added.
! 27 October 2008. Version 1.4.0.  Optional argument n_threads added.

! This is based on the code:
! -- A Fully Portable High Performance Minimal --
! -- Storage Hybrid Format Cholesky Algorithm  --
!    B.S. Andersen - DMU, Denmark;
!    J.A. Gunnels & F.G. Gustavson - IBM, USA;
!    J.K. Reid - RAL, UK;
!    J. Wasniewski - DTU, Denmark
!    December 10, 2004.

module hsl_ma54_double

  implicit none
  EXTERNAL dgemm, dtrsm, dsyrk, dcopy, dtpsv, dgemv

  integer,parameter,private :: wp = kind(1.0d0) ! Precision parameter.
  integer,parameter,private :: long = selected_int_kind(18) ! Long integer.
  real(wp),parameter,private :: unity = 1.0_wp
  integer(long),parameter,private :: one = 1_long

  interface ma54_to_block
      module procedure ma54_to_block_double
  end interface

  interface ma54_from_block
      module procedure ma54_from_block_double
  end interface

  interface ma54_factor
      module procedure ma54_factor_double
  end interface

  interface ma54_solve1
      module procedure ma54_solve1_double
  end interface

  interface ma54_solve2
      module procedure ma54_solve2_double
  end interface

  interface ma54_forward1
      module procedure ma54_forward1_double
  end interface

  interface ma54_forward2
      module procedure ma54_forward2_double
  end interface

  interface ma54_back1
      module procedure ma54_back1_double
  end interface

  interface ma54_back2
      module procedure ma54_back2_double
  end interface

  interface ma54_diag
      module procedure ma54_diag_double
  end interface

contains

subroutine ma54_to_block_double(n,p,nb,ap,buf,info,partial)
! Rearrange a packed triangular matrix to blocked hybrid format
   integer, intent (in) :: n ! Specifies the matrix order
   integer, intent (in) :: p ! Column p is be at the end of a block
   integer, intent (in) :: nb ! Block size for the blocked hybrid format.
   real (wp), intent (inout) :: ap((n*(n+one))/2)
!                Holds the matrix, which is rearranged to the new format.
   real (wp) :: buf(nb*int(n,long)) ! Work array.
   integer, intent (out) :: info ! set to unity of these values:
!                               0 Successful solution.
!                             < 0 Failure
   logical, optional, intent(in) :: partial ! Ignored.

!  Local variables
   integer :: jb ! Length of block
   integer :: j1 ! First column of block
   integer :: i, j ! Row, column index within the block column
   integer(long) :: ijap, ijbuf ! Positions of entry in ap and buf
   integer(long)  :: ac ! Position in ap of start of block column
   integer :: l ! Length of row or column
   integer :: lb ! Length of block column

   intrinsic min

   info = 0
   if (p > n)  info = -3
   if (n < 0)  info = -1
   if (p < 0)  info = -2
   if (nb < 1) info = -5
   if (info/=0 .or. n==0) return

   ijap = 1
! Main loop over block columns
   do j1 = 1,p,nb
      jb = min(nb,p-j1+1)
      lb = n - j1 + 1
      ac = ijap
! Copy to the buffer
      ijbuf = 1
      do j = 1, jb
         l = lb - j + 1
         call dcopy(l,ap(ijap),1,buf(ijbuf),1)
         ijbuf = ijbuf + lb + 1
         ijap = ijap + l
      end do

! Copy array back from buffer
      ijap = ac
      do i = 1, lb
         l = min(i,jb)
         call dcopy(l,buf(i),lb,ap(ijap),1)
        ijap = ijap + l
      end do
   end do

end subroutine ma54_to_block_double

subroutine ma54_factor_double(n,p,nb,ap,buf,info,n_threads)
! Cholesky part factorization of a symmetric positive-definite matrix held
! in packed blocked hybrid format.

!$     use omp_lib

  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Factorize to column p
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (inout) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix
!              A in packed blocked hybrid format. The ma54_factor of A is
!              contained in packed blocked hybrid format on exit.
  real (wp) :: buf(nb*int(n,long)) ! Work array.
  integer, intent (out) :: info ! set to unity of these values:
!                       0 Successful solution.
!                    /= 0 Failure:
!                     > 0. The leading minor of order info is not positive
!                       definite, and the factorization could not be completed.
   integer, intent (inout), optional :: n_threads ! Controls the number of
!              threads used in OpenMP parallel regions. Its value may be changed
!              by another thread during the execution of this subroutine.
!              0  Execute on omp_get_num_threads() threads.
!              1  Execute on a single thread without reference to OpenMP.
!             >1  Execute on n_threads threads.
!             Absence is treated as equivalent to presence with the value 0.

!  Local variables
    integer i ! Temporary variable
    integer ib ! Length of block row
    integer(long) ii ! Position in ap of diagonal entry
    integer j ! First column of current block column
    integer jb ! Number of columns in current block column
    integer(long) jbt ! Size of packed triangle of order jb
    integer(long) jj ! Position in ap of start of current diagonal block
    integer(long) jk ! Position in ap of start of current block
    integer kb ! Number of columns in earlier block column
    integer k ! First column of current block
    integer l ! Length of column
    integer lb ! Length of block column
    integer(long) kk ! Position in buf of diagonal entry
    integer nt  ! Number of threads

    real (wp) ddot
    intrinsic min

    info = 0
    if (p > n)  info = -3
    if (n < 0)  info = -1
    if (p < 0)  info = -2
    if (nb <= 0) info = -5
    if (info/=0 .or. n==0) return

    jj = 0
! Loop over leading block columns
    do j = 1, p, nb
      jb = min(nb,p-j+1)
      jbt = (jb*(jb+one))/2
      ii = jj
      kk = 1
! Move diagonal block to buf in full upper format
      do i = 1, jb
! Transfer i-th row of ap to i-th col of buf, if n > nb
        call dcopy(i,ap(ii),1,buf(kk),1)
        ii = ii + i
        kk = kk + jb
      end do

! Choose the number of threads for the next parallel loop
      nt = 1
!$    nt = omp_get_max_threads()
      if( nt>1 ) then
         if ( present(n_threads) ) then
            if( n_threads>0 ) nt = n_threads
         end if
      end if
!$    if (nt>0) call omp_set_num_threads(nt)

! Apply the previous elimination operations to the block column
!$omp parallel default(none) &
!$omp private(i, k, jk, kb, ib) &
!$omp shared(ap, j, jb, n, nb, jj, jbt, p, buf, info, lb)
! Note that the implied barrier at end of single avoids write-write
! data races in the syrk updates, allowing a NOWAIT on the gemm loop
      jk = 0
      do k = 1, min(j-1,p), nb ! Loop over block columns
        kb = min(nb,p-k+1)
        jk = jk + (j-k)*int(kb,long) - (kb*(kb-one))/2
! Update the diagonal block, held in buf
!$omp single
        call dsyrk('U','T',jb,kb,-unity,ap(jk),kb,unity,buf,jb)
!$omp end single
!$omp do
        do i = j+jb, n, nb
           ib = min(n+1-i, nb)
           call dgemm('T','N',jb,ib,kb,-unity,ap(jk),kb, &
              ap(jk+kb*int(i-j,long)),kb,unity, &
              ap(jj+jbt+jb*int(i-j-jb)),jb)
        end do
!$omp end do nowait
        jk = jk + (n+one-j)*kb
      end do
! Cholesky factorization of diagonal block
!$omp single
        call dkcf(jb,buf,jb,i)
!       call dpotrf('u',jb,buf,jb,i)
        if (i>0) info = i + j - 1
!$omp end single
! Important: implied barrier at end of single flushes info
        if(info.eq.0) then
! Calculate the sub-diagonal blocks of L in the current block column
!$omp do
          do i = j+jb, n, nb
            ib = min(n+1-i, nb)
            call dtrsm('L','U','T','N',jb,ib,unity,buf,jb, &
              ap(jj+jbt+jb*int(i-j-jb)),jb)
          end do
!$omp end do
        end if
!$omp end parallel
      if(info.ne.0) return

      ii = jj
      kk = 1
! Move  diagonal block back from buf to ap
      do i = 1, jb
        call dcopy(i,buf(kk),1,ap(ii),1)
        ii = ii + i
        kk = kk + jb
      end do

      jj = jj + jb*(n-j+one) - (jb*(jb-one))/2
    end do

! Loop over trailing block columns
    do j = p+1, n, nb
      jb = min(nb,n-j+1)
      jbt = (jb*(jb+one))/2
      ii = jj
      kk = 1
! Move block column to buf in rectangular format
      lb = n - j + 1
      do i = 1, jb
        l = lb - i + 1
        call dcopy(l,ap(ii),1,buf(kk),1)
        ii = ii + l
        kk = kk + lb + 1
      end do

! Choose the number of threads for the next parallel loop
      nt = 1
!$    nt = omp_get_max_threads()
      if( nt>1 ) then
         if ( present(n_threads) ) then
            if( n_threads>0 ) nt = n_threads
         end if
      end if
!$    if (nt>0) call omp_set_num_threads(nt)

! Apply the previous elimination operations to the block column
!$omp parallel default(none) &
!$omp private(i, k, jk, kb, ib) &
!$omp shared(ap, j, jb, n, nb, jj, jbt, p, buf, info, lb)
! Note that the implied barrier at end of single avoids write-write
! data races in the syrk updates, allowing a NOWAIT on the gemm loop
      jk = 0
      do k = 1, min(j-1,p), nb ! Loop over block columns
        kb = min(nb,p-k+1)
        jk = jk + (j-k)*int(kb,long) - (kb*(kb-one))/2
! Update the diagonal block, held in buf
!$omp single
        call dsyrk('L','T',jb,kb,-unity,ap(jk),kb,unity,buf,lb)
!$omp end single
! Update the off-diagonal part of the block column
!$omp do
        do i = j+jb, n, nb
          ib = min(n+1-i, nb)
          call dgemm('T','N',ib,jb,kb,-unity, &
            ap(jk+kb*int(i-j,long)),kb,ap(jk),kb, unity,buf(i-j+1),lb)
        end do
!$omp end do nowait

        jk = jk + (n+one-j)*kb
      end do
!$omp end parallel
      if(info.ne.0) return

      ii = jj
      kk = 1
! Move block column back from buf to ap
      do i = 1, jb
         l = lb - i + 1
         call dcopy(l,buf(kk),1,ap(ii),1)
         ii = ii + l
         kk = kk + lb + 1
      end do

      jj = jj + jb*(n-j+one) - (jb*(jb-one))/2
    end do

end subroutine ma54_factor_double


subroutine ma54_from_block_double(n,p,nb,ap,buf,info,partial)
! Rearrange blocked dhybrid packed matrix to packed triangular format
   integer, intent (in) :: n ! Specifies the matrix order
   integer, intent (in) :: p ! Column p is at the end of a block
   integer, intent (in) :: nb ! Block size  for the blocked hybrid format.
   real (wp), intent (inout) :: ap((n*(n+one))/2) ! Holds the matrix, which
!                       is rearranged to the packed triangular format.
   real (wp) :: buf(nb*int(n,long)) ! Work array.
   integer, intent (out) :: info ! set to unity of these values:
!                           0 Successful solution.
!                         < 0 Failure:
   logical, optional, intent(in) :: partial ! Ignored

   integer :: jb ! Column length of block
   integer :: j1 ! First column of block
   integer :: i, j ! Row, column index within the block
   integer(long) :: ijap, ijbuf ! Positions of ap(i1+i-1,j1+j-1), buf(i,j)
   integer(long) :: ac ! Position in ap of start of block
   integer :: l ! Length of row or column
   integer :: lb ! Length of block column
   intrinsic min

   info = 0
   if (p > n)  info = -3
   if (n < 0)  info = -1
   if (p < 0)  info = -2
   if (nb < 1) info = -5
   if (info/=0 .or. n==0) return

   ijap = 1
! Main loop over block columns
   do j1 = 1, p, nb
      jb = min(nb,p-j1+1)
      lb = n - j1 + 1

! Copy to buffer
      ac = ijap
      do i = 1, lb
         l = min(i,jb)
         call dcopy(l,ap(ijap),1,buf(i),lb)
         ijap = ijap + l
      end do

! Copy from the buffer
      ijap = ac
      ijbuf = 1
      do j = 1, jb
         l = lb - j + 1
         call dcopy(l,buf(ijbuf),1,ap(ijap),1)
         ijbuf = ijbuf + lb + 1
         ijap = ijap + l
      end do
   end do

end subroutine ma54_from_block_double


subroutine ma54_forward2_double(n,p,nbo,nrhs,ap,b,ldb,mbo,buf,info)
! Partial forward substitution for unity or more sets of equations, given a
!                partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nbo ! Block size for the blocked hybrid format.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  real (wp), intent (inout) :: b(ldb,*) ! Holds the right-hand sides on
!                              entry and is overwritten by the solution.
  integer, intent (in) :: mbo ! Block size for the right-hand sides.
  real (wp) :: buf(int(p,long)*nbo + int(n,long)*mbo)
                              ! Work array (needed if nrhs >=4)
  integer, intent (out) :: info ! set to unity of these values:
!                         0 Successful solution.
!                       < 0 Failure

! Local variables
  integer(long) :: bufb ! Position in buf that holds a buffer for b
  integer :: i ! Row index
  integer :: ib ! Number of columns in current block column of b
  integer(long) :: ibp ! Current position in bufb
  integer(long) :: ii ! Position in ap of diagonal entry
  integer(long) :: ij ! Position in ap
  integer :: j ! Column index
  integer :: jb ! Number of columns in current block column of ap
  integer(long) :: jbp ! Current position in bufb
  integer(long) :: jbd ! nb*int(kb,long), increment for jbp
  integer(long) :: jd ! size of current trapezoid / increment for jj
  integer(long) :: jj ! Position in ap of start of current diagonal block
  integer :: k ! First column of current block
  integer :: kb ! ! Number of columns in current block column of b
  integer(long) :: kk ! Position in buf of diagonal entry
  integer(long) :: kks ! Current position in buf
  integer :: mb ! block size used for the blocked rhs:min(m,mbo)
  integer :: nb ! block size used for the blocked hybrid format: min(n,nbo)
  integer(long) :: nb2 ! nb*int(nb,long)
  integer(long) :: nbt ! Size of packed triangle of order nb
  integer(long) :: nb1 ! Size of packed triangle of order nb-1

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nrhs<0) info = -4
  if (nbo<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mbo<1) info = -7
  nb = min(n,nbo)
  mb = min(mbo,nrhs)
  if (info/=0 .or. n==0) return

  if (nrhs<4) then
    do j = 1, nrhs
      call ma54_forward1_double(n,p,nbo,ap,b(1,j),info)
    end do
    return
  end if

  bufb = p*int(nb,long)

  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  nb1 = nbt - nb

!  Load the diagonal blocks (triangles) of ap into buf.
  kks = 1
  jd = n*int(nb,long) - nb1
  jj = 0
  do j = 1, p, nb
    jb = min(p-j+1,nb)
    ii = jj
    kk = kks
    do i = 1, jb
      call dcopy(i,ap(ii),1,buf(kk),1)
      ii = ii + i
      kk = kk + jb
    end do
    jj = jj + jd ! -> next triangle in ap
    jd = jd - nb2 ! next delta for jj
    kks = kks + nb2 ! -> next triangle of ap in buf
  end do ! do j

! Main loop over the block columns of b
  do k = 1, nrhs, mb
    kks = 1
    jd = n*int(nb,long) - nb1
    jj = 0
    kb = min(nrhs-k+1,mb)

! Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
    jbp = 1
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do

! Solve U'*Y=B
    jbd = nb*int(kb,long) ! delta for jbp
    jbp = 1
    do j = 1, p, nb
      jb = min(p-j+1,nb)
      call dtrsm('L','U','T','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      ibp = jbp + jb*int(kb,long) ! -> starting block of B
      ij = jj + (jb*(jb+one))/2 ! -> starting block of AP
      do i = j + nb, p, nb
        ib = min(p-i+1,nb)
        call dgemm('T','N',ib,kb,jb,-unity,ap(ij),jb,buf(bufb+jbp),jb,unity, &
          buf(bufb+ibp),ib)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      do i = p+1, n, nb
        ib = min(n-i+1,nb)
        call dgemm('T','N',ib,kb,jb,-unity,ap(ij),jb,buf(bufb+jbp),jb,unity, &
          buf(bufb+ibp),ib)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      jj = jj + jd
      jd = jd - nb2
      jbp = jbp + jbd ! -> next block of B
      kks = kks + nb2 ! -> next triangle of ap
    end do ! do j

!       Copy solution to b(1:n,k:k+kb+1)
    jbp = 1 ! -> to buf
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
  end do ! do k

end subroutine ma54_forward2_double

subroutine ma54_forward1_double(n,p,nb,ap,b,info)
! Partial forward substitution for unity set of equations, given a
!          partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) ! b(1:n) holds the right-hand
!                   sides on entry and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!          0 Successful solution.
!        < 0 The value of argument -info is not valid.

! Local variables
  integer j ! First column of current block column
  integer jb ! Number of columns in current block column
  integer(long) jbt ! Size of packed triangle of order jb
  integer(long) jd ! size of current trapezoid / increment for jj
  integer(long) jj ! Position in ap of start of current diagonal block
  integer(long) nb2 ! nb*int(nb,long)
  integer(long) nbt ! Size of packed triangle of order nb

! .. Executable Statements ..
  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nb<1) info = -5

  if (info/=0 .or. n==0)  return
  jj = 0
  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  jd = n*int(nb,long) - nbt + nb

!     Solve U'*Y=B
  do j = 1, p - nb, nb
    call dtpsv('U','T','N',nb,ap(jj),b(j),1)
    call dgemv('T',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j),1,unity,b(j+nb),1)
    jj = jj + jd
    jd = jd - nb2
  end do
  jb = p - j + 1
  jbt = (jb*(jb+one))/2
  call dtpsv('U','T','N',jb,ap(jj),b(j),1)
  if (n>=j+jb) then
    call dgemv('T',jb,n+1-j-jb,-unity,ap(jj+jbt),jb,b(j),1,unity,b(j+jb),1)
  end if
end subroutine ma54_forward1_double


subroutine ma54_back2_double(n,p,nb,nrhs,ap,b,ldb,mb,buf,info)
! Partial back substitution for unity or more sets of equations, given a
!              partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  real (wp), intent (inout) :: b(ldb,*) ! Holds the right-hand sides on
!                             entry and is overwritten by the solution.
  integer, intent (in) :: mb ! Block size for the right-hand sides.
  real (wp) :: buf(p*int(nb,long) + n*int(mb,long)) ! Work array
  integer, intent (out) :: info ! set to unity of these values:
!                 0 Successful solution.
!               < 0 Failure

! Local variables
  integer(long) :: bufb
  integer :: i ! Row index
  integer :: ib ! Number of columns in current block column of b
  integer(long) :: ibp ! Current position in bufb
  integer(long) :: ii ! Position in ap of diagonal entry
  integer(long) :: ij ! Position in ap
  integer :: j ! Column index
  integer :: jb ! Number of columns in current block column of ap
  integer(long) :: jbp ! Current position in bufb
  integer(long) :: jbd ! nb*int(kb,long), increment for jbp
  integer(long) :: jd ! size of current trapezoid / increment for jj
  integer(long) :: jj ! Position in ap of start of current diagonal block
  integer :: k ! First column of current block
  integer :: kb ! ! Number of columns in current block column of b
  integer(long) :: kk ! Position in buf of diagonal entry
  integer(long) :: kks ! Current position in buf
  integer(long) :: nb2 ! nb*int(nb,long)
  integer(long) :: nbt ! Size of packed triangle of order nb
  integer(long) :: nb1 ! Size of packed triangle of order nb-1

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nrhs<0) info = -4
  if (nb<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mb<1) info = -7
  if (info/=0 .or. n==0) return

  if (nrhs<4) then
    do j = 1, nrhs
      call ma54_back1_double(n,p,nb,ap,b(1,j),info)
    end do
    return
  end if

  bufb = p*int(nb,long)

  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  nb1 = nbt - nb

!  Load the diagonal blocks (triangles) of ap into buf.
  kks = 1
  jd = n*int(nb,long) - nb1
  jj = 0
  do j = 1, p, nb
    jb = min(p-j+1,nb)
    ii = jj
    kk = kks
    do i = 1, jb
      call dcopy(i,ap(ii),1,buf(kk),1)
      ii = ii + i
      kk = kk + jb
    end do
    jj = jj + jd ! -> next triangle in ap
    jd = jd - nb2 ! next delta for jj
    kks = kks + nb2 ! -> next triangle of ap in buf
  end do ! do j

!     Main loop over the block columns of b
  do k = 1, nrhs, mb
    kb = min(nrhs-k+1,mb)

! Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
    jbp = 1
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
        jbp = jbp + ib
      end do
    end do

! Solve U*X=Y
    j = ((p-1)/nb)*int(nb,long)
    jj = (j*(n+n-j+one))/2
    jbp = 1+j*int(kb,long)
    jbd = nb*int(kb,long) ! delta for jbp
    jd = (nb*(2*(n-j)+nb+one))/2
    kks = 1+j*int(nb,long)
    do j = j+1, 1, -nb
      jb = min(p-j+1,nb)
      ij = jj + (jb*(jb+one))/2 ! -> starting block of AP
      ibp = jbp + jb*int(kb,long) ! -> starting update block of B
      do i = j + nb, p, nb
        ib = min(p-i+1,nb)
        call dgemm('N','N',jb,kb,ib,-unity,ap(ij),jb,buf(bufb+ibp),ib,unity, &
          buf(bufb+jbp),jb)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long) ! -> next block of AP
      end do ! do i
      do i = p+1, n, nb
        ib = min(n-i+1,nb)
        call dgemm('N','N',jb,kb,ib,-unity,ap(ij),jb,buf(bufb+ibp),ib,unity, &
          buf(bufb+jbp),jb)
        ibp = ibp + ib*int(kb,long) ! -> next block of B
        ij = ij + ib*int(jb,long)  ! -> next block of AP
      end do ! do i
      call dtrsm('L','U','N','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      jj = jj - jd ! -> next triangle of AP
      jd = jd + nb2 ! next delta for jj
      jbp = jbp - jbd ! -> current block of B
      kks = kks - nb2 ! -> current triangle of ap
    end do ! do j

! Copy solution to b(1:n,k:k+kb+1)
    jbp = 1 ! -> to buf
    do i = 1, p, nb
      ib = min(p-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
    do i = p+1, n, nb
      ib = min(n-i+1,nb)
      do j = k, k + kb - 1
        call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
        jbp = jbp + ib
      end do
    end do
  end do ! do k

end subroutine ma54_back2_double

subroutine ma54_back1_double(n,p,nb,ap,b,info)
! Partial back substitution for unity set of equations, given a
!             partial Cholesky factorization in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(p*(n*2-p+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) ! b(1:n) holds the right-hand
! sides on entry and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!                   0 Successful solution.
!                 < 0 The value of argument -info is not valid.

  integer i ! Temporary variable
  integer j ! First column of current block column
  integer jb ! Number of columns in current block column
  integer(long) jbt ! Size of packed triangle of order jb
  integer(long) jd ! size of current trapezoid / increment for jj
  integer(long) jj ! Position in ap of start of current diagonal block
  integer(long) nb2 ! nb*int(nb,long)
  integer(long) nbt ! Size of packed triangle of order nb

  info = 0
  if (p>n) info = -3
  if (n<0) info = -1
  if (p<0) info = -2
  if (nb<1) info = -5
  if (info/=0 .or. n==0) return

  jj = 0
  nb2 = nb*int(nb,long)
  nbt = (nb2+nb)/2
  jd = n*int(nb,long) - nbt + nb
  do j = 1, p - nb, nb
    jj = jj + jd
    jd = jd - nb2
  end do
  jb = p - j + 1
  jbt = (jb*(jb+one))/2
  if (n>=j+jb) then
    call dgemv('N',jb,n+1-j-jb,-unity,ap(jj+jbt),jb,b(j+jb),1,unity,b(j),1)
  end if
  call dtpsv('U','N','N',jb,ap(jj),b(j),1)
  i = j - nb
  do j = i, 1, -nb
    jd = jd + nb2
    jj = jj - jd
    call dgemv('N',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j+nb),1,unity,b(j),1)
    call dtpsv('U','N','N',nb,ap(jj),b(j),1)
  end do

end subroutine ma54_back1_double

subroutine ma54_solve1_double(n,nb,ap,b,info)
! Solve unity or more sets of equations, given the Cholesky factorization
!                  of its matrix in blocked hybrid format.
! For the character argument, the case is insignificant.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix.
  real (wp), intent (inout) :: b(n) !  holds the right-hand side on entry
!            and is overwritten by the solution.
  integer, intent (out) :: info ! set to unity of these values:
!                0 Successful solution.
!              < 0 Failure

    integer i ! Temporary variable
    integer j ! First column of current block column
    integer jb ! Number of columns in current block column
    integer(long) jd ! size of current trapezoid / increment for jj
    integer(long) jj ! Position in ap of start of current diagonal block
    integer(long) nb2 ! nb*int(nb,long)
    integer(long) nbt ! Size of packed triangle of order nb

  info = 0
  if (n<0) info = -1
  if (nb<1) info = -5

  if (info/=0 .or. n==0) return

    jj = 0
    nb2 = nb*int(nb,long)
    nbt = (nb2+nb)/2
    jd = n*int(nb,long) - nbt + nb

!     Solve U'*Y=B
    do j = 1, n - nb, nb
      call dtpsv('U','T','N',nb,ap(jj),b(j),1)
      call dgemv('T',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j),1,unity,b(j+nb),1)
      jj = jj + jd
      jd = jd - nb2
    end do
    jb = n - j + 1
    call dtpsv('U','T','N',jb,ap(jj),b(j),1)

!     Solve U*X=Y
    call dtpsv('U','N','N',jb,ap(jj),b(j),1)
    i = j - nb
    do j = i, 1, -nb
      jd = jd + nb2
      jj = jj - jd
      call dgemv('N',nb,n+1-j-nb,-unity,ap(jj+nbt),nb,b(j+nb),1,unity,b(j),1)
      call dtpsv('U','N','N',nb,ap(jj),b(j),1)
    end do
 end subroutine ma54_solve1_double

subroutine ma54_solve2_double(n,nb,nrhs,ap,b,ldb,mb,buf,info)
! Solve unity or more sets of equations, given the Cholesky factorization
!              of its matrix in blocked hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap(0:(n*(n+one))/2-1) ! Holds the matrix.
  integer, intent (in) :: ldb ! The first dimension of the array b.
  integer, intent (in) :: nrhs ! The number of right-hand sides.
  real (wp), intent (inout) :: b(ldb,*) ! b(1:n,1:nrhs) holds the
!             right-hand sides on entry and is overwritten by the solution.
  integer, intent (in) :: mb ! Block size for the right-hand sides.
  real (wp) :: buf(*) ! work array
  integer, intent (out) :: info ! set to unity of these values:
!                                0 Successful solution.
!                              < 0 Failure

    integer(long) :: bufb
    integer :: i ! Row index
    integer :: ib ! Number of columns in current block column of b
    integer(long) :: ibp ! Current position in bufb
    integer(long) :: ii ! Position in ap of diagonal entry
    integer(long) :: ij ! Position in ap
    integer :: j ! Column index
    integer :: jb ! Number of columns in current block column of ap
    integer(long) :: jbp ! Current position in bufb
    integer(long) :: jbd ! nb*int(kb,long), increment for jbp
    integer(long) :: jd ! size of current trapezoid / increment for jj
    integer(long) :: jj ! Position in ap of start of current diagonal block
    integer :: k ! First column of current block
    integer :: kb ! ! Number of columns in current block column of b
    integer(long) :: kk ! Position in buf of diagonal entry
    integer(long) :: kks ! Current position in buf
    integer(long) :: nb2 ! nb*int(nb,long)
    integer(long) :: nbt ! Size of packed triangle of order nb
    integer(long) :: nb1 ! Size of packed triangle of order nb-1
    intrinsic min

  info = 0
  if (n<0) info = -1
  if (nrhs<0) info = -4
  if (nb<1) info = -5
  if (ldb<max(1,n)) info = -6
  if (mb<1) info = -7

  if (info/=0 .or. n==0 .or. nrhs==0) return

! If there are few rhs, solve each separately
    if (nrhs<4) then
      do j = 1, nrhs
        call ma54_solve1_double(n,nb,ap,b(1,j),info)
      end do
      return
    end if

    bufb = n*int(nb,long)

    nb2 = nb*int(nb,long)
    nbt = (nb2+nb)/2
    nb1 = nbt - nb

!     Load the diagonal blocks (triangles) of ap into buf.
    kks = 1 ! note that kks is not restored to 1 here.
    jd = n*int(nb,long) - nb1 ! jd is not restored to this value here.
    jj = 0 ! is not restored to zero here.
    do j = 1, n, nb
      jb = min(n-j+1,nb)
      ii = jj
      kk = kks
      do i = 1, jb
        call dcopy(i,ap(ii),1,buf(kk),1)
        ii = ii + i
        kk = kk + jb
      end do
      jj = jj + jd ! -> next triangle in ap
      jd = jd - nb2 ! next delta for jj
      kks = kks + nb2 ! -> next triangle of ap in buf
    end do ! do j

!     Main loop over the block columns of b
    kks = 1 ! is changed but restored to 1
    jd = n*int(nb,long) - nb1 ! is changed but restored to this value
    jj = 0 ! is changed but restored to 0
    do k = 1, nrhs, mb
      kb = min(nrhs-k+1,mb)

!       Copy kb RHS from  b(1:n,k:k+kb-1) to buffer bufb
      jbp = 1 ! -> to bufb
      do i = 1, n, nb
        ib = min(n-i+1,nb)
        do j = k, k + kb - 1
          call dcopy(ib,b(i,j),1,buf(bufb+jbp),1)
          jbp = jbp + ib
        end do
      end do

!      Solve U'*Y=B
      jbd = nb*int(kb,long) ! delta for jbp
      jbp = 1 ! is changed but restored to 1
      do j = 1, n - nb, nb
        call dtrsm('L','U','T','N',nb,kb,unity,buf(kks),nb,buf(bufb+jbp),nb)
        ibp = jbp + jbd ! -> starting block of B
        ij = jj + nbt ! -> starting block of AP
        do i = j + nb, n, nb
          jb = min(n-i+1,nb)
          call dgemm('T','N',jb,kb,nb,-unity,ap(ij),nb,buf(bufb+jbp),nb,unity, &
            buf(bufb+ibp),jb)
          ibp = ibp + jbd ! -> next block of B
          ij = ij + nb2 ! -> next block of AP
        end do ! do i
        jj = jj + jd
        jd = jd - nb2
        jbp = jbp + jbd ! -> next block of B
        kks = kks + nb2 ! -> next triangle of ap
      end do ! do j
      jb = n - j + 1
      call dtrsm('L','U','T','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)

!       Solve U*X=Y
      call dtrsm('L','U','N','N',jb,kb,unity,buf(kks),jb,buf(bufb+jbp),jb)
      do j = j - nb, 1, -nb
        jd = jd + nb2
        jj = jj - jd ! -> next triangle of AP
        ij = jj + nbt ! -> starting block of AP
        ibp = jbp ! -> starting update block of B
        jbp = jbp - jbd ! -> current block of B
        kks = kks - nb2 ! -> current triangle of ap
        do i = j + nb, n, nb
          jb = min(n-i+1,nb)
          call dgemm('N','N',nb,kb,jb,-unity,ap(ij),nb,buf(bufb+ibp),jb,unity, &
            buf(bufb+jbp),nb)
          ibp = ibp + jbd ! -> next block of B
          ij = ij + nb2 ! -> next block of AP
        end do ! do i
        call dtrsm('L','U','N','N',nb,kb,unity,buf(kks),nb,buf(bufb+jbp),nb)
      end do ! do j

!       Copy solution to b(1:n,k:k+kb+1)
      jbp = 1 ! -> to buf
      do i = 1, n, nb
        ib = min(n-i+1,nb)
!         copy buf(is:is+ib+k-1) to B(i:i+ib-1,0:k-1); K <= NB
        do j = k, k + kb - 1
          call dcopy(ib,buf(bufb+jbp),1,b(i,j),1)
          jbp = jbp + ib
        end do
      end do
    end do ! do k

  end subroutine ma54_solve2_double

 subroutine ma54_diag_double(n,p,nb,ap,diag)
! Extract the diagonal of a partial Cholesky factorization in blocked
! hybrid format.
  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Number of columns
  integer, intent (in) :: nb ! Block size for the blocked hybrid format.
  real (wp), intent (in) :: ap((n*(n+one))/2) ! Holds the matrix.
  real (wp), intent (out) :: diag(p) ! Diagonal of factorization.

  integer js ! First column of current block column
  integer j  ! Column index within the block
  integer(long) k  ! Index within ap
  integer l  ! Length of rectangular part of the block column

  k = 0
  l = n-nb
  do js = 1, p, nb
     do j = 1, min(nb,p-js+1)
       k = k + j
       diag(js+j-1) = ap(k)
     end do
     k = k + l*int(nb,long)
     l = l - nb
  end do
end subroutine ma54_diag_double


subroutine dkcf(n,a,lda,info)
! Kernel subroutine for Cholesky factorization
  integer, intent (in) :: n ! Specifies the matrix order.
  integer, intent (in) :: lda ! Leading extent of array a.
! The inequality lda >= n must hold.
! For efficiency, the value n for lda is preferable.
  real (wp), intent (inout) :: a(0:lda-1,0:n-1) ! The upper-triangular part,
! a(i,j), i<=j<=n, must be set to hold the upper-triangular part
! of the matrix and is overwritten by the Cholesky ma54_factor.
! For efficient execution, a(:,1:n) should fit into
! level-1 cache.
  integer, intent (out) :: info ! set to unity of these values:
! 0 Successful solution.
! i>0. Pivot i is not positive.
! .. Locals ..
  real (wp), parameter :: zero = 0.0_wp, unity = 1.0_wp
  integer :: i ! Temporary variable
  integer :: j ! Temporary variable
  integer :: k ! Temporary variable
  integer :: nr ! Temporary variable
  real (wp) :: rd0 ! scalar to hold reciprocal of diagonal A(K,K) for 0<=K<N
  real (wp) :: t00, t01, t02, t03 ! regs to hold a(j,i:i+3)
  real (wp) :: ai0, ai1, ai2, ai3 ! regs to hold a(k,i:i+3)
!     or A(K:K+3,J) or A(K:K+3,J+1)

!     this code uses 8 FP regs



!     This is a 1 x 4 blocking kernel
!     It is ONLY optimized for the case MOD(N,4) = 0

  info = 0
  nr = mod(n,4)

!     Main loop: ma54_factor A(0:n-nr-1,0:n-nr-1) = U'*U

  do j = 0, n - nr - 4, 4

!       Process row J = A(J,J:N-1)

!       Process A(J  ,J  )

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j)
      ai1 = a(k+1,j)
      ai2 = a(k+2,j)
      ai3 = a(k+3,j)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 2
    ai0 = sqrt(ai0)
    a(j,j) = ai0
    rd0 = unity/ai0


!       process 1 by 3 rectangle A(J,J+1:J+3)

    t00 = a(j,j+1)
    t01 = a(j,j+2)
    t02 = a(j,j+3)
    do k = 0, j - 1
      ai0 = a(k,j+1)
      ai1 = a(k,j+2)
      ai2 = a(k,j+3)
      ai0 = ai0*a(k,j)
      ai1 = ai1*a(k,j)
      ai2 = ai2*a(k,j)
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
    end do

!       scale and store A(J,J+1:J+3)

    t00 = t00*rd0
    t01 = t01*rd0
    t02 = t02*rd0
    a(j,j+1) = t00
    a(j,j+2) = t01
    a(j,j+3) = t02

!       Process the remainder of rows J = A(J,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j,i)
      t01 = a(j,i+1)
      t02 = a(j,i+2)
      t03 = a(j,i+3)
      do k = 0, j - 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j)
        ai1 = ai1*a(k,j)
        ai2 = ai2*a(k,j)
        ai3 = ai3*a(k,j)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j,i) = t00
      a(j,i+1) = t01
      a(j,i+2) = t02
      a(j,i+3) = t03
    end do

!       Process row J+1 = A(J+1,J+1:N-1)

!       Process A(J+1,J+1)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+1)
      ai1 = a(k+1,j+1)
      ai2 = a(k+2,j+1)
      ai3 = a(k+3,j+1)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+1)
    ai0 = ai0*ai0
    t00 = t00 - ai0

    ai0 = a(j+1,j+1)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 3
    ai0 = sqrt(ai0)
    a(j+1,j+1) = ai0
    rd0 = unity/ai0


!       process 1 by 2 rectangle A(J+1,J+2:J+3)

    t00 = a(j+1,j+2)
    t01 = a(j+1,j+3)
    do k = 0, j
      ai0 = a(k,j+2)
      ai1 = a(k,j+3)
      ai0 = ai0*a(k,j+1)
      ai1 = ai1*a(k,j+1)
      t00 = t00 - ai0
      t01 = t01 - ai1
    end do

!       scale and store A(J+1,J+2:J+3)

    t00 = t00*rd0
    t01 = t01*rd0
    a(j+1,j+2) = t00
    a(j+1,j+3) = t01

!       Process the remainder of row J+1 = A(J+1,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+1,i)
      t01 = a(j+1,i+1)
      t02 = a(j+1,i+2)
      t03 = a(j+1,i+3)
      do k = 0, j
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+1)
        ai1 = ai1*a(k,j+1)
        ai2 = ai2*a(k,j+1)
        ai3 = ai3*a(k,j+1)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+1,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+1,i) = t00
      a(j+1,i+1) = t01
      a(j+1,i+2) = t02
      a(j+1,i+3) = t03
    end do

!       Process row J+2 = A(J+2,J+2:N-1)

!       Process A(J+2,J+2)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+2)
      ai1 = a(k+1,j+2)
      ai2 = a(k+2,j+2)
      ai3 = a(k+3,j+2)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+2)
    ai1 = a(j+1,j+2)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    t00 = t00 - ai0
    t01 = t01 - ai1

    ai0 = a(j+2,j+2)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 4
    ai0 = sqrt(ai0)
    a(j+2,j+2) = ai0
    rd0 = unity/ai0


!       process 1 by 1 rectangle A(J+2,J+3)

    t00 = a(j+2,j+3)
    do k = 0, j + 1
      ai0 = a(k,j+3)
      ai0 = ai0*a(k,j+2)
      t00 = t00 - ai0
    end do

!       scale and store A(J+2,J+3)

    t00 = t00*rd0
    a(j+2,j+3) = t00

!       Process the remainder of row J+2 = A(J+2,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+2,i)
      t01 = a(j+2,i+1)
      t02 = a(j+2,i+2)
      t03 = a(j+2,i+3)
      do k = 0, j + 1
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+2)
        ai1 = ai1*a(k,j+2)
        ai2 = ai2*a(k,j+2)
        ai3 = ai3*a(k,j+2)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+2,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+2,i) = t00
      a(j+2,i+1) = t01
      a(j+2,i+2) = t02
      a(j+2,i+3) = t03
    end do

!       Process row J+3 = A(J+3,J+3:N-1)

!       Process A(J+3,J+3)

    t00 = zero
    t01 = zero
    t02 = zero
    t03 = zero

    do k = 0, j - 1, 4
      ai0 = a(k,j+3)
      ai1 = a(k+1,j+3)
      ai2 = a(k+2,j+3)
      ai3 = a(k+3,j+3)
      ai0 = ai0*ai0
      ai1 = ai1*ai1
      ai2 = ai2*ai2
      ai3 = ai3*ai3
      t00 = t00 - ai0
      t01 = t01 - ai1
      t02 = t02 - ai2
      t03 = t03 - ai3
    end do
    ai0 = a(j,j+3)
    ai1 = a(j+1,j+3)
    ai2 = a(j+2,j+3)
    ai0 = ai0*ai0
    ai1 = ai1*ai1
    ai2 = ai2*ai2
    t00 = t00 - ai0
    t01 = t01 - ai1
    t02 = t02 - ai2

    ai0 = a(j+3,j+3)
    t00 = t00 + t01
    t02 = t02 + t03
    t00 = t00 + t02
    ai0 = ai0 + t00

    if (ai0<=zero) go to 5
    ai0 = sqrt(ai0)
    a(j+3,j+3) = ai0
    rd0 = unity/ai0


!       Process the remainder of row J+3 = A(J+3,J+4:N-1)
!       using 1 by 4 rectangles

    do i = j + 4, n - 4, 4
      t00 = a(j+3,i)
      t01 = a(j+3,i+1)
      t02 = a(j+3,i+2)
      t03 = a(j+3,i+3)
      do k = 0, j + 2
        ai0 = a(k,i)
        ai1 = a(k,i+1)
        ai2 = a(k,i+2)
        ai3 = a(k,i+3)
        ai0 = ai0*a(k,j+3)
        ai1 = ai1*a(k,j+3)
        ai2 = ai2*a(k,j+3)
        ai3 = ai3*a(k,j+3)
        t00 = t00 - ai0
        t01 = t01 - ai1
        t02 = t02 - ai2
        t03 = t03 - ai3
      end do

!         scale and store A(J+3,I+0:I+3)

      t00 = t00*rd0
      t01 = t01*rd0
      t02 = t02*rd0
      t03 = t03*rd0
      a(j+3,i) = t00
      a(j+3,i+1) = t01
      a(j+3,i+2) = t02
      a(j+3,i+3) = t03
    end do
  end do
  if (nr==0) return

!     here 0 < NR < 4
!     scale last NR cols of A

  call dtrsm('L','U','T','N',n-nr,nr,unity,a,lda,a(0,n-nr),lda)

!     update triangle A(N-NR:N-1,N-NR,N-1) with A(0:N-NR-1,N-NR:N-1)

  call dsyrk('U','T',nr,n-nr,-unity,a(0,n-nr),lda,unity,a(n-nr,n-nr),lda)

!     Cholesky ma54_factor triangle A(N-NR:N-1,N-NR,N-1)

  j = n - nr
  if (nr==3) then ! ma54_factor row J and update A(J+1:J+2,J+1:J+2)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = unity/t00
    ai0 = a(j,j+1)
    ai1 = a(j,j+2)
    ai0 = ai0*rd0
    ai1 = ai1*rd0
    a(j,j+1) = ai0
    a(j,j+2) = ai1
    t00 = a(j+1,j+1)
    t01 = a(j+1,j+2)
    t02 = a(j+2,j+2)
    ai2 = ai0*ai0
    ai3 = ai0*ai1
    ai1 = ai1*ai1
    t00 = t00 - ai2
    t01 = t01 - ai3
    t02 = t02 - ai1
    a(j+1,j+1) = t00
    a(j+1,j+2) = t01
    a(j+2,j+2) = t02
    j = j + 1
  end if
  if (nr>=2) then ! ma54_factor row J and update A(J+1,J+1)
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
    rd0 = unity/t00
    ai0 = a(j,j+1)
    ai0 = ai0*rd0
    a(j,j+1) = ai0
    t00 = a(j+1,j+1)
    ai1 = ai0*ai0
    t00 = t00 - ai1
    a(j+1,j+1) = t00
    j = j + 1
  end if
  if (nr>=1) then ! A(J,J) = SQRT( A(J,J) )
    t00 = a(j,j)
    if (t00<=zero) go to 2
    t00 = sqrt(t00)
    a(j,j) = t00
  end if

  return

2 info = j + 1
  return
3 info = j + 2
  return
4 info = j + 3
  return
5 info = j + 4
  return
end subroutine dkcf

end module hsl_ma54_double
! COPYRIGHT (c) 2007 Council for the Central Laboratory
!               of the Research Councils
!
! Version 6.2.0
! For version history see ChangeLog
!
! To change precision:
!    Change dgemm, dcopy, dswap, dgemv, dtpsv, daxpy
!    Change _double, kind(1.0d0)

module hsl_ma64_double

  implicit none
  EXTERNAL dgemm, dcopy, dswap, dgemv, dtpsv, daxpy
  private
  public ma64_factor, ma64_solveL1, ma64_solveD1, &
         ma64_solveDLT1, ma64_solveLT1, &
         ma64_solveL2, ma64_solveD2, ma64_solveDLT2, &
         ma64_solveLT2, ma64_control, ma64_info
  integer, parameter :: wp = kind(1.0d0) ! Precision parameter.
  integer,parameter:: long = selected_int_kind(18) ! Long integer.
  real (wp), parameter :: unity = 1.0_wp

  interface ma64_factor
      module procedure ma64_factor_double
  end interface

  interface ma64_solveL1
     module procedure ma64_solveL1_double
  end interface

  interface ma64_solveD1
     module procedure ma64_solveD1_double
  end interface

  interface ma64_solveDLT1
     module procedure ma64_solveDLT1_double
  end interface

  interface ma64_solveLT1
     module procedure ma64_solveLT1_double
  end interface

  interface ma64_solveL2
     module procedure ma64_solveL2_double
  end interface

  interface ma64_solveD2
     module procedure ma64_solveD2_double
  end interface

  interface ma64_solveDLT2
     module procedure ma64_solveDLT2_double
  end interface

  interface ma64_solveLT2
     module procedure ma64_solveLT2_double
  end interface

  type ma64_control
     integer :: p_thresh=32 ! If p<=p_thresh, execute on a single thread.
     real (wp) :: small=1e-20_wp ! Diagonal entries of D of modulus less
!                      than this are replaced by zero.
     real (wp) :: static=0.0_wp ! If static>0 and p stable pivots are not found,
!                      the 1x1 pivot that is nearest to satisfying the test is
!                      chosen and info%num_nothresh is incremented by one.
!                      If its absolute value is less than cntl%static, it is
!                      replaced by the nearer of -static and static and
!                      info%num_perturbed is incremented by one.
     logical :: twos=.false. ! If true, the signs of perm indicate recommended
!                      2x2 pivots.
     real (wp) :: u=0.1_wp ! Initial relative pivot tolerance. Values greater
!                      than 0.5 are treated as 0.5 and values less than zero are
!                      treated as zero.
     real (wp) :: umin=1.0_wp ! Minimum relative pivot tolerance. Values greater
!                      than u are treated as u and values less than zero are
!                      treated as zero. If p stable pivots have not been not
!                      found and the candidate pivot with greatest relative
!                      pivot tolerance has tolerance v>=umin, this is accepted
!                      as a pivot and u is reset to v. If p=n> and both u and
!                      umin are greater than 0.5, umin is treated as having the
!                      value 0.5.
end type ma64_control

  type ma64_info
     real (wp) :: detlog=0 ! log of the absolute value of the determinant
!                          of D or zero if the determinant is zero.
     integer :: detsign=0 !  the sign of the determinant of D or zero
!                         if the determinant of D is zero.
!                         Not used in the complex symmetric case.

     integer :: flag=0 ! Zero after successful entry, negative after error.
     integer :: num_neg=0 ! Number of negative eigenvalues of D.
!                           Not used in the complex symmetric case.
     integer :: num_nothresh=0 !  the number diagonal entries of D that
!                      were chosen as 1x1 pivots without satisfying the
!                      relative pivot threshold.
     integer :: num_perturbed=0 ! Number diagonal entries of D that
!                      were perturbed to cntl%static or -cntl%static.
     integer :: num_zero=0 ! Number of zero eigenvalues of D.
     integer :: num_2x2=0 ! Number of 2x2 blocks in D.
     real (wp) :: usmall ! Set to -1 if num_perturbed > 0. Otherwise,
!              if q = p, it holds the smallest relative pivot value of the
!                  chosen pivots, or
!              if q < p and a positive value of cntl%umin would have led to a
!                  greater value of q, it holds the largest such value;
!                  otherwise, it holds zero.
     real (wp) :: u=0 ! Set to the final value of the relative pivot
!                 threshold u.
  end type ma64_info


contains

subroutine ma64_factor_double &
             (n,p,nb,nbi,a,la,cntl,q,ll,perm,d,buf,info,s,n_threads)
! Partial factorization of a symmetric indefinite matrix held in lower
! packed format to the form
!                 P A P' = ( L   ) ( D   ) (L' M')
!                          ( M  I) (   S ) (   I )
! where P is a permutation matrix that involves only the first p rows and
! columns, L is unit lower triangular of size q and D is block diagonal of
! order q with blocks of size 1 and 2.

!$     use omp_lib

  integer, intent (in) :: n ! Specifies the matrix order
  integer, intent (in) :: p ! Pivots to be chosen in rows and cols 1:p.
  integer, intent (in) :: nb ! Block size
  integer, intent (in) :: nbi ! Inner block size
  integer(long), intent (in) :: la ! Size of a. la>= min (n*n,(n*(n+nb+1))/2).
  real (wp), intent (inout) :: a(la) ! On entry, a(la+1-ln:la) holds the
!              matrix in lower packed format, where ln=(n*(n+1))/2. On exit,
!              a(1:ll) holds ( L ) as a sequence of block columns, each of which
!                            ( M )
!              consists of a triangular matrix packed by columns followed by a
!              rectangular matrix packed by columns. There are nb columns in
!              each block column except the last, which may have fewer. Also,
!              a(la+1-ls:la) holds the matrix S in lower packed format,
!              where ls=((n-q)*(n-q+1))/2.
  type(ma64_control), intent (in) :: cntl ! Controls the actions.
  integer, intent (out) :: q ! Order of D.
  integer(long), intent (out) :: ll ! Size of (L M)' in block-column form,
!              that is, (q*(n+n-q+1))/2.
  integer, intent (inout) :: perm(p) ! The input value is ignored if
!              cntl%twos is false; otherwise, each sequence
!              perm(i)<0, perm(i+1)<0, ... perm(i+k)<0 must be of even length
!              and flags recommended 2x2 pivots.
!              On return, for i = 1, 2, ..., p, perm(i) is set to the index
!              of the row of A that was permuted to row i.
  real (wp), intent (out) :: d(2*p) ! d(1:2*q) is set to hold the inverse of
!              D, except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2). d(2*q) is set to zero.
  real (wp) :: buf(nb*n+n) ! Work array for the current part of -LD.
  type(ma64_info), intent (inout) :: info
  integer, intent (in), optional :: s ! Columns 1:s should be searched last for
!              pivots.
  integer, intent (inout), optional, target :: n_threads ! Controls the number
!              of threads used in OpenMP parallel regions. Its value may be
!              changed by another thread during the execution of this
!              subroutine.
!              0  Execute on omp_get_num_threads() threads.
!             >1  Execute on n_threads threads.
!             Absence is treated as equivalent to presence with the value 0.

!     .. Local Scalars ..
    real (wp) amax  ! Entry of largest absolute value to left of diagonal in
!                     row m.
    real (wp) amax2 ! Second largest absolute value to left of diagonal in row m
    real (wp) amaxb ! Largest absolute value below diagonal in column m
    real (wp) amaxt ! Largest absolute value in column t.
    integer   deti   ! Determinant held as detr*radix**deti
    real (wp) detpiv ! Determinant of candidate pivot is detpiv/detscale.
    real (wp) detpiv0 ! First term in calculation of detpiv.
    real (wp) detpiv1 ! Second term in calculation of detpiv.
    real (wp) detpiv2 ! Value of detpiv for best candidate 2x2 pivot.
    real (wp) detr   ! Scaled determinant
    real (wp) detscale ! Inverse of the largest entry in the candidate pivot.
    real (wp) detscale2 ! Value of detscale for best candidate 2x2 pivot.
    integer i ! Row index
    integer i1 ! First row of block
    integer j ! Column index
    integer j1 ! First column of block
    integer j2 ! Last column of block
    integer js ! Last column available for swapping
    integer(long) k  ! Position in a
    integer(long) kdiag  ! Position in a of diagonal np+1
    integer(long) kkj ! Position of diagonal of column j
    integer(long) kkt ! Position of diagonal of column t
    integer(long) kkr1 ! Position of diagonal of column r+1
    integer(long) kkq !  Start of the pivotal block column
    integer(long) kkq1 ! Position of diagonal of column q+1
    integer(long) kkm ! Position of diagonal of column m
    integer(long) kkn ! Position of diagonal of column np
    integer(long) kkmb ! Position of diagonal of column mbest
    integer(long) kkm1 ! Position of diagonal of column m-1
    integer(long) kq1 ! Position in buf of diagonal of column q+1
    integer lj ! Height of block containing column j
    integer lm ! Height of block containing column m
    integer lq ! Height of the pivot block
    integer m  ! Column searched most recently
    integer m1  ! First column of the block in which m appears
    integer mbest ! Column with best relative pivot value for a 1x1 pivot
    integer mbest1, mbest2 ! Pair with best relative pivot value for a 2x2 pivot
    integer mdummy  ! Loop execution count
    integer mlast ! Last column of the block in which m appears
    integer np ! Columns 1 to np are rearranged to block column form.
!                np = p or np = n.
    integer nt  ! Number of threads
    integer pivsiz ! Size of the chosen pivot, 0 if none chosen, or -1
!                    if column is essentially zero
    integer q1 ! First column of the block containing column q+1
    integer q2 ! Last column of the inner block containing column q+1
    integer r ! Number of pivot operations applied to columns q+1:p
    integer rm ! Number of pivot operations applied to the (outer) block
               ! containing column m
    real (wp) rmax ! Largest entry in column m
    real (wp) rmax1 ! Largest entry in column m-1 outwith rows m,m-1 if (m,m-1)
                    ! is a recommended 2x2 pivot. Otherwise, -1.
    real (wp) rmax2 ! Largest entry in column m outwith rows m,m-1.
    integer t ! Candidate 2x2 pivot is in columns t and m
    real (wp) ubest1 ! Relative pivot value of best candidate 1x1 pivot
    real (wp) ubest2 ! Relative pivot value of best candidate 2x2 pivot
    real (wp) urel ! Relative pivot value
    real (wp) urel1 ! Relative pivot value in column m-1 if (m,m-1) is a
!                    recommended 2x2 pivot.
    real (wp) u ! Relative pivot threshold
    real (wp) umin ! Minimum relative pivot threshold

    info%flag = 0
    if (n < 0) then
       info%flag = -1
    else if (p < 0) then
       info%flag = -2
    else if (p > n) then
       info%flag = -3
    else if (nbi <= 1) then
       info%flag = -4
    else if ( la < min(n*int(n,long),(n*(n+nb+1_long))/2) ) then
       info%flag = -7
    else if (cntl%static < cntl%small .and. cntl%static/=0.0_wp) then
       info%flag = -10
    else if (mod(nb,nbi) /= 0) then
       info%flag = -12
    end if
    info%detlog = 0.0_wp
    info%num_neg = 0
    info%num_nothresh = 0
    info%num_perturbed = 0
    info%num_zero = 0
    info%num_2x2 = 0
    u = min(max(cntl%u,0.0_wp),1.0_wp)
    info%u = u
    umin = min(max(cntl%umin,0.0_wp),u)
    if (p==n) umin = min(umin,0.5_wp)
    info%usmall = 1.0_wp
    q = 0
    if (info%flag/=0 .or. p==0) return
    deti = 0
    detr = 1.0_wp
    nt = 1
!$  if (p>cntl%p_thresh) nt = omp_get_max_threads()
    np = n
    if( nt==1 .and. p<=nb ) np = p
    lq = n
    kkq1 = 1
    kkq = 1
    m = p
! m is updated at the start of the main loop so initializing it to p causes
! it to have the value 1 during the first execution of that loop.
    js = p
    q1 = 1
    q2 = min(nbi,p)
    rmax1 = -1.0_wp

! For j = 1,p: set perm(j) = -j for the first part of a 2x2 pivot; otherwise,
! set perm(j) = j
    if (cntl%twos) then
      js = 0
      do j = 1, p
        if (perm(j) >= 0) then
          perm(j) = j
        else
! First part of a 2x2 pivot
          perm(j) = -j
          if (j==p) then
             info%flag = -11
             return
          else if (perm(j+1) >= 0) then
             info%flag = -11
             return
          else
             perm(j+1) = j+1
          end if
        end if
      end do
    else
      do j = 1, p
        perm(j) = j
      end do
    end if

! Rearrange first np columns to block column form
    call to_block

  if (present(s)) then
     if (s>0 .and. s<p) then
        js = js - s
        i = min(s,p-s)
! Ignore recommended 2x2 pivots that will be split
        if (perm(i)<0) perm(i) = i
        if (perm(p-i)<0) perm(p-i) = p-i
! Make first i columns be last
        do j = 1,i
           call swap_cols(j,p+1-j)
        end do
! Ensure that the first index of a 2x2 pair is labelled
        do j = 2,i
           if (perm(j)<0) then
               perm(j) = -perm(j)
               perm(j-1) = -perm(j-1)
           end if
        end do
        do j = p+2-i,p
           if (perm(j)<0) then
               perm(j) = -perm(j)
               perm(j-1) = -perm(j-1)
           end if
        end do
     end if
  end if

  pivot: do
! Perform a pivotal operation
    pivsiz = 0
    ubest1 = 0.0_wp
    ubest2 = 0.0_wp
    sweep: do mdummy = 1, p-q ! Look for a pivotal column or pair of columns

! Update m and the scalars associated with column m
      m = m+1
      if (m>p) then
! Go back to column q+1
         m = q+1
         m1 = q1
         kkm = kkq1
         lm = lq
         r = q
         rm = q
      else if (m<m1+nb) then
! Within the current block column
         kkm = kkm + lm + 1
      else
! At the the start of another block column
         lm = lm - nb
         m1 = m1 + nb
         rm = r
         kkm = kkm + lm + 1
      end if

! If there are more than nbi active columns and there are enough columns
! not yet considered, update all candidate columns and perform swaps.
      if (m-q>nbi) then
        if (js-m > m-q) then
           mlast = min(p,m1+nb-1)
!     Inner update: of columns (m:mlast) for pivots rm+1:q (BLAS3)
           if (m<=mlast .and. rm<q) call update (m,mlast,rm+1,q)
!     Outer update: of columns (m1+nb:p) for pivots r+1:q (BLAS3)
           if (m1+nb<=p .and. r<q) call update (m1+nb,p,r+1,q)
           r = q
           rm = r
           do j = q+1,m
              call swap_cols(j,js)
              js = js -1
           end do
        end if
      end if

! Update column m
      kkr1 = kkq1 + (rm-q)*(lq+1_long)
      k = kkr1+m-rm-1
      if(q>rm) then
         call dgemv('NoTrans',n-m+1,q-rm,unity,a(k),lq, &
                    buf(n*(rm+1_long-q1)+m),n,unity,a(kkm),1)
      end if


! Find largest and second largest entry to left of diagonal in row m.
       j = q + 1
       lj = lq
       k = kkq1 + m - j - lj
       amax = 0.0_wp
       amax2 = 0.0_wp
       t = 0
       if (j<m) then
          t = j
          kkt = kkq1
       end if
       do j1 = q1,m-1,nb
          do j = j, min(j1+nb-1,m-1)
             k = k + lj
             if (abs(a(k))>abs(amax)) then
                t = j
                amax2 = abs(amax)
                amax = a(k)
                kkt = k - (m-j)
             else
                amax2 = max(abs(a(k)),amax2)
             end if
          end do
          lj = lj - nb
       end do

! Now calculate largest entry below the diagonal of column m.
       amaxb = 0.0_wp
       do i = m+1,n
          amaxb = max(abs(a(kkm+i-m)),amaxb)
       end do

! Now calculate largest entry in the whole of column m and make sure that is
! is neither small nor infinity.
       rmax = max(abs(amax),abs(a(kkm)),amaxb)
       if (rmax<=cntl%small) then
! All entries of the column are small
           a(kkm) = 0.0_wp
           info%num_zero = info%num_zero + 1
           pivsiz = -1
           perm(m) = abs(perm(m))
           rmax1 = -1.0_wp
           exit sweep
       else if (rmax > huge(1.0_wp)) then
! There is an infinity in the column
           info%flag = -13
           return
       end if

! Calculate the relative pivot value and see if it is the best so far
       if (abs(a(kkm))>cntl%small) then
          urel = abs(a(kkm))/rmax
       else
          urel = 0.0_wp
       end if
       if (urel >= ubest1) then
          ubest1 = urel
          mbest = m
          kkmb = kkm
       end if

! If (m,m+1) recommended as 2x2 pivot store needed data and cycle so that
! column m+1 is updated
       if (perm(m)<0) then
          rmax1 = abs(amax)
          do i = m+2,n
             rmax1 = max(abs(a(kkm+i-m)),rmax1)
          end do
          urel1 = urel
          kkm1 = kkm
          perm(m) = abs(perm(m))
          cycle sweep
       end if

! If (m-1,m) recommended as 2x2 pivot, look for pivot in columns m-1,m
       lookm1m: if (rmax1>=0.0_wp) then

! Find largest entry in column m outwith rows m-1 and m
          rmax2 = amaxb
          j = q + 1
          lj = lq
          k = kkq1 + m - j - lj
          do j1 = q1,m-2,nb
             do j = j, min(j1+nb-1,m-2)
                k = k + lj
                rmax2 = max(abs(a(k)),rmax2)
             end do
             lj = lj - nb
          end do
          k = kkm1 + 1

! If both diagonal entries >= cntl%small or the off-diagonal entry
! >= cntl%small, consider 2x2 pivot
          test2x2: if ( min(abs(a(kkm)),abs(a(kkm1))).ge.cntl%small .or. &
              abs(a(k)).ge.cntl%small ) then
             detscale = 1.0_wp/max( abs(a(kkm)), abs(a(kkm1)), abs(a(k)) )
             detpiv1 =  (abs(a(k))*detscale)*abs(a(k))
             detpiv0 = (a(kkm)*detscale)*a(kkm1)
             detpiv = detpiv0 - detpiv1
! Make sure that the 2x2 pivot is not singular and that there is little
! cancellation in calculating detpiv. Bearing in mind the way detscale
! is calculated, if the largest entry of the matrix is chosen as pivot,
! the one entry of the reduced matrix has absolute value abs(detpiv).
             if (abs(detpiv)> &
                 max(cntl%small,abs(detpiv0)/2,abs(detpiv1)/2)) then
! OK as 2x2 pivot if all entries in the rest of the columns are small
                if(max(rmax2,rmax1)<=cntl%small) then
                   pivsiz = 2
                   rmax1 = -1.0_wp
                   t = m-1
                   exit sweep
                end if

! Calculate the relative pivot value (as 2x2 pivot)
                urel = abs(detpiv)/max(                         &
                    abs(a(kkm)*detscale)*rmax1+abs(a(k)*detscale)*rmax2,   &
                    abs(a(kkm1)*detscale)*rmax2+abs(a(k)*detscale)*(rmax1) )
                rmax1 = -1.0_wp

! OK as 2x2 pivot if relative pivot value is big enough
                if (urel>u) then
                   pivsiz = 2
                   t = m-1
                   info%usmall = min(urel,info%usmall)
                   exit sweep
                end if

! If this has the best relative pivot value so far, record this
                if(urel>ubest2)then
                   ubest2 = urel
                   detpiv2 = detpiv
                   detscale2 = detscale
                   mbest1 = m
                   mbest2 = m-1
                end if

             end if

          end if test2x2

! If m-1 OK as 1x1 pivot, accept this
          rmax1 = -1.0_wp
          if (urel1>u) then
             if (abs(a(kkm1))>cntl%small) then
                pivsiz = 1
                info%usmall = min(urel1,info%usmall)
                call swap_cols(m-1,m)
                exit sweep
             end if
          end if

       end if lookm1m

! If there is a candidate 2x2 pivot, try it.
! Look for a 1x1 pivot only if the 2x2 pivot is unacceptable.
       tgt0: if (t>0) then
       if ( min(abs(a(kkm)),abs(a(kkt))).ge.cntl%small .or. &
              abs(amax).ge.cntl%small ) then
! Store value of the largest entry in whole of column m outwith rows m and t
          rmax2 = max(amax2,amaxb)

! Find largest entry in column t outwith rows m and t
          amaxt = 0.0_wp
          j = q + 1
          lj = lq
          k = kkq1 + t - j - lj
          do j1 = q1,t-1,nb
            do j = j, min(j1+nb-1,t-1)
               k = k + lj
               amaxt = max(abs(a(k)),amaxt)
            end do
            lj = lj - nb
          end do
          if (j1/=t) lj = lj + nb
          k = k + lj
          do i = t+1,m-1
            amaxt = max(abs(a(k+i-t)),amaxt)
          end do
          do i = m+1,n
            amaxt = max(abs(a(k+i-t)),amaxt)
          end do

           detscale = 1.0_wp/max( abs(a(kkm)), abs(a(kkt)), abs(amax) )
          detpiv1 =  (amax*detscale)*amax
          detpiv0 = a(kkm)*detscale*a(kkt)
          detpiv = detpiv0 - detpiv1
! Make sure that the 2x2 pivot is not singular and that there is little
! cancellation in calculating detpiv. Bearing in mind the way detscale
! is calculated, if the largest entry of the matrix is chosen as pivot,
! the one entry of the reduced matrix has absolute value abs(detpiv).
          left2x2:if (abs(detpiv)> &
              max(cntl%small,abs(detpiv0)/2,abs(detpiv1)/2)) then

! OK as 2x2 pivot if all entries in the rest of the columns are small
            if(max(rmax2,amaxt)<=cntl%small) then
                pivsiz = 2
                exit sweep
            end if

! Calculate the relative pivot value (as 2x2 pivot)
            urel = abs(detpiv)/max( &
               abs(a(kkm)*detscale)*amaxt+(abs(amax)*detscale)*rmax2, &
               abs(a(kkt)*detscale)*rmax2+(abs(amax)*detscale)*amaxt )

! OK as 2x2 pivot if relative pivot value is big enough
            if (urel>u) then
              pivsiz = 2
               info%usmall = min(urel,info%usmall)
               exit sweep
            end if

! If this has the best relative pivot value so far, record this
            if(urel>ubest2)then
               ubest2 = urel
               detpiv2 = detpiv
               detscale2 = detscale
               mbest1 = m
               mbest2 = t
            end if
         end if left2x2
      end if
      end if tgt0

! If 2x2 pivot rejected or only one column left, take best 1x1 pivot if is OK.
      if (t>0 .or. m==p) then
         if (ubest1>u) then
            pivsiz = 1
            info%usmall = min(ubest1,info%usmall)
            if (mbest/=m) call swap_cols(mbest,m)
            exit sweep
         end if
      end if

    end do sweep

! No pivot found in search of all available columns
    pivsiz0: if (pivsiz==0) then
! Since all the columns have been updated, revise m to q+1
        m = q+1
        kkm = kkq1
        m1 = q1
        lm = lq
        rm = q
        r = q
! Perform relaxed pivoting if the best pivot is good enough
        if (max(ubest1,ubest2)>=umin) then
           if (ubest1>=ubest2) then
! Accept 1x1 pivot
              u = min(u,ubest1)
              info%usmall = min(ubest1,info%usmall)
              pivsiz = 1
              if (mbest/=m) call swap_cols(m,mbest)
           else
! Accept 2x2 pivot
              u = min(u,ubest2)
              info%usmall = min(ubest2,info%usmall)
              pivsiz = 2
! Revise m to q+1
              m = q+2
              kkm = kkq1+lq+1
              if (m==m1+nb) then
                 m1 = m1+nb
                 lm = lm - nb
              end if
              detpiv = detpiv2
              detscale = detscale2
              t = min(mbest1,mbest2)
              if (t/=q+1) call swap_cols(q+1,t)
              t = max(mbest1,mbest2)
              if (t/=q+2) call swap_cols(q+2,t)
              t = q+1
           end if

! Perform static pivoting if this has been requested
        else if (cntl%static>0) then
           if (mbest/=m) call swap_cols(m,mbest)
           pivsiz = 1
           info%num_nothresh = info%num_nothresh + 1
           if (abs(a(kkm)) < cntl%static) then
               a(kkm) = sign(cntl%static,real(a(kkm),wp))
               info%num_perturbed = info%num_perturbed + 1
               info%usmall = -1
           else
               info%usmall = min(ubest1,info%usmall)
           end if
        end if
    end if pivsiz0

    pivsizes: if (pivsiz==1) then
! Swap columns if m not q+1
      if (q+1/=m) call swap_cols(q+1,m)
! We store D**-1.
      d(2*q+1) = 1.0_wp/a(kkq1)
      d(2*q+2) = 0.0_wp
! Update info
!      info%detlog = info%detlog + log(abs(a(kkq1)))
      deti = deti + exponent(detr) + exponent(abs(a(kkq1)))
      detr = fraction(detr)*fraction(abs(a(kkq1)))
      if (a(kkq1)<0.0_wp) info%num_neg = info%num_neg + 1
! Store L in A and -LD or -conjg(LD) in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = -a(kkq1)
      a(kkq1) = 1.0_wp
      do i = 1, n-q-1
         buf(kq1+i) = -a(kkq1+i)
         a(kkq1+i) = d(2*q+1)*a(kkq1+i)
      end do
      j1 = q1
      lj = lq
! Update columns q+2 to m
      j = q+2
      kkj = kkq1
      do j1 = j1,m,nb
         do j = j, min(j1+nb-1,m)
            kkj = kkj + lj + 1
            call daxpy(n-j+1,buf(kq1+j-q-1),a(kkq1+j-q-1),1,a(kkj),1)
         end do
         lj = lj - nb
      end do
! Update q and kkq1
      kkq1 = kkq1 + lq + 1
      if(q+2==q1+nb)kkq1 = kkq1 - nb
      q = q+1

    else if (pivsiz==2) then pivsizes
! Swap columns unless t==q+1 and m==q+2
      if (q+2/=m) then
         if (q+1/=t) call swap_cols(q+1,t)
         call swap_cols(q+2,m)
      end if
! We store D**-1
      k = kkq1 + lq + 1
      if(q+2==q1+nb)k = k - nb
      d(2*q+1) = (a(k)*detscale)/detpiv
      d(2*q+3) = (a(kkq1)*detscale)/detpiv
      d(2*q+2) = -(a(kkq1+1)*detscale)/detpiv
      d(2*q+4) = 0.0_wp

! Update info
      info%num_2x2 = info%num_2x2 + 1
!     info%detlog = info%detlog + log(abs(detpiv)) - log(detscale)
      deti = deti + exponent(detr) - exponent(detscale)
      detr = fraction(detr)*abs(detpiv)/fraction(detscale)
      if (detpiv<0.0_wp) then
         info%num_neg = info%num_neg + 1
      else if (a(kkq1)+a(k)<0.0_wp) then
         info%num_neg = info%num_neg + 2
      end if
! Store L in A and -LD or -conjg(LD) in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = -a(kkq1)
      buf(kq1+1) = -a(kkq1+1)
      buf(kq1+n) = -a(kkq1+1)
      buf(kq1+n+1) = -a(k)
      a(kkq1) = 1.0_wp
      a(k) = 1.0_wp
      a(kkq1+1) = 0.0_wp
      do i = 2, n-q-1
         buf(kq1+i  ) = -a(kkq1+i)
         buf(kq1+i+n) = -a(k+i-1)
         a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
         a(k+i-1) = d(2*q+3)*a(k+i-1) - d(2*q+2)*buf(kq1+i)
      end do
! Update columns q+3 to m
      j1 = q1
      lj = lq
      if(q+3>=q1+nb) then
! Pivot is at end of block
         lj = lq - nb
         j1 = q1 + nb
      end if
      kkj = k + lj + 1
      do j1 = j1,m,nb
         j = max(q+3,j1)
         j2 = min(j1+nb-1,m)
         if(j2>=j) then
            call dgemm('n','t',n-j+1,j2-j+1,2,unity, &
                   a(kkq1+j-q-1),n-((q+1)/nb)*nb,buf(kq1+j-q-1),n, &
                   unity,a(kkj),lj)
         end if
         kkj = kkj + (lj+1_long)*(j2-j+1) - nb
         lj = lj - nb
      end do
! Update q and kkq1
      kkq1 = k + lq + 1
      if (q+3>=q1+nb) kkq1 = kkq1 - nb ! Pivot is at end of block
      q = q + 2

    else if (pivsiz==-1) then pivsizes
! Handle a row that is zero
      if (q+1/=m) call swap_cols(q+1,m)
      d(2*q+1) = 0.0_wp
      d(2*q+2) = 0.0_wp
! Store L in A and -LD in buf
      kq1 = q+1+(q+1_long-q1)*n
      buf(kq1) = 0.0_wp
      a(kkq1) = 1.0_wp
      do i = 1, n-q-1
         buf(kq1+i) = 0.0_wp
         a(kkq1+i) = 0.0_wp
      end do
      kkq1 = kkq1 + lq + 1
      if(q+2==q1+nb)kkq1 = kkq1 - nb
      q = q+1
    end if pivsizes

    j = m1+nb
    if (q>=q2) then
!     Inner block is complete
!     See if the inner update can be done with the outer update
      if( q>=q1+nb-1 .and. r==rm ) then
         j = m+1
      else
!     Update columns (m+1:min(p,m1+nb-1)) for pivots rm+1:q2 (BLAS3)
        mlast = min(p,m1+nb-1)
        if (m<mlast) call update (m+1,mlast,rm+1,q2)
      end if
      rm = q2
      q2 = min (p,q2+nbi)
    end if

    if (q>=q1+nb-1) then
!     Outer block is complete
!     Update columns (j:p) for pivots r+1:q1+nb-1 (BLAS3)
      if (j<=p) call update (j,p,r+1,q1+nb-1)
      r = q1+nb-1
      rm = r

      if (p<n) then
!     Update columns (p+1:n) for pivots q1:q1+nb-1
         if(np<n) then
            call update2 (q1+nb-1)
         else
            call update (p+1,n,q1,q1+nb-1)
         end if
      end if

      q1 = q1 + nb
      if (q1==q) then
!     If a 2x2 pivot spans the blocks, move final column of buf forward
         call dcopy(n-q+1,buf(q+n*int(nb,long)),1,buf(q),1)
      end if
      kkq = kkq + int(nb,long)*lq
      lq = lq - nb

    end if

    if (q == p .or. pivsiz==0) exit pivot
  end do pivot

  if (q>=q1 .and. p<n) then
!  Update columns (p+1:n) for pivots q1:q (BLAS3)
     if(np<n) then
        call update2 (q)
     else
        call update (p+1,n,q1,q)
     end if
  end if

! Update info
  if(q<p .and. info%num_nothresh==0) &
            info%usmall = min(info%usmall, max(ubest1,ubest2))
  if(info%num_zero/=0) then
     info%detsign = 0
     info%detlog = 0.0_wp
  else
     info%detsign = -1
     if (mod(info%num_neg,2)==0)info%detsign = 1
     info%detlog = log(detr)+deti*log(1.0d0*radix(1.0d0))
  end if
  info%u = u

! Rearrage columns 1 to np
  call from_block




contains

    subroutine update(jl,jr,pl,pr)
!     Update columns jl:jr for pivots pl:pr when jr<=np.
       integer jl,jr,pl,pr
! Local variables
       integer ij ! Index for looping over blocks
       integer i2 ! Last row of block
       integer j  ! Column index
       integer nc  ! Number of block columns to be updated
       integer nr  ! Number of block rows to be updated
       integer nn  ! Number of blocks to be updated
       integer j0  ! Last column of block ahead of jl
       integer p0  ! First column of pivot block

       j0 = nb*((jl-1)/nb)
       p0 = nb*((pl-1)/nb)+1
       nr = 1 + (n-j0-1)/nb
       nc = 1 + (jr-j0-1)/nb
       nn = (nc*(2*nr-nc+1))/2

! Decide on the number of threads to use
       nt = 1
       if (nn>1 .and. p>cntl%p_thresh .and. np==n) then
!$        nt = omp_get_max_threads()
          if( nt>1 ) then
             if ( present(n_threads) ) then
                if( n_threads>0 ) nt = n_threads
             end if
          end if
       end if

       if(nt==1) then  ! No parallelization
          do j1 = j0+1, jr, nb
             j = max(jl,j1)
             j2 = min(j1+nb-1,jr)
             if(n.ge.j .and. j2.ge.j) then
                kkj = j1-1
                kkj = 1 + ((kkj)*(n+n-kkj+nb))/2
                if(j.ne.j1) kkj = kkj + (j-j1)*(n-j0+1)
                call dgemm &
                    ('n','t',n-j+1,j2-j+1,pr-pl+1,unity, &
                    a(kkq+j-q1+lq*(pl-p0)),lq,buf(j+n*(pl-p0)), &
                    n,unity,a(kkj),n-j1+1)
            endif
          end do

       else ! Parallelize by blocks
!$        if (nt>0) call omp_set_num_threads(nt)
!$OMP     PARALLEL DO &
!$OMP     DEFAULT(SHARED) PRIVATE(I, J, I1, J1, I2, J2, KKJ) &
!$OMP     SCHEDULE(STATIC)
          do ij = 1, nn
! Find the block indices
             j = nr + 1.4999d0 - sqrt( (nr+0.5d0)**2 -ij*2 )
             i = ij - ((j-1)*(nr+nr-j))/2
! Find the leading and trailing row and column
             i1 = j0+1+(i-1)*nb
             j1 = j0+1+(j-1)*nb
             i2 = min(i1+nb-1,n)
             j2 = min(j1+nb-1,jr)
! Find the first row and column to be updated
             j = max(jl,j1)
             i = max(jl,i1)
! Find position of diagonal of column j
             kkj = j1-1
             kkj = 1 + ((kkj)*(n+n-kkj+nb))/2
             if(j/=j1) kkj = kkj + (j-j1)*(n-j0+1)
             call dgemm &
                 ('n','t',i2-i+1,j2-j+1,pr-pl+1,unity,a(kkq+i-q1+lq*(pl-p0)) &
                        ,lq,buf(j+n*(pl-p0)),n,unity,a(kkj+i-j),n-j1+1)
          end do
!OMP      END PARALLEL DO
       end if
    end subroutine update

    subroutine update2(pr)
!     Update columns p+1:n for pivots q1:pr within the current block column
!     when np<n. Columns p+1:n are held in packed format.
      integer pr
! Local variables
       integer j1  ! First column of block column
       integer j2  ! Last column of block column
       integer(long) k ! Column start in packed form
       integer(long) kfree ! Start of free space in a
       integer(long) l ! Column start in unpacked form
       integer p0  ! First column of pivot block

       p0 = nb*((q1-1)/nb)+1
       kkj = kdiag
       kfree = kkn + n+1-np
       do j1 = p+1, n, nb/2
          j2 = min(j1+nb/2-1,n)
! Copy block of columns to rectangular format
          k = kkj
          l = kfree
          do j = j1, j2
             call dcopy(n+1-j,a(k),1,a(l),1)
             k = k + n+1-j
             l = l + n+2-j1
          end do
! Perform the update
          call dgemm &
              ('n','t',n+1-j1,j-j1,pr-q1+1,unity,a(kkq+j1-q1+lq*(q1-p0)), &
                     lq,buf(j1+n*(q1-p0)),n,unity,a(kfree),n+1-j1)
! Copy block of columns back
          k = kkj
          l = kfree
          do j = j1, j2
             call dcopy(n+1-j,a(l),1,a(k),1)
             k = k + n+1-j
             l = l + n+2-j1
          end do
          kkj = k
       end do
    end subroutine update2

    subroutine to_block
! Rearrange the first np columns of a packed triangular matrix to block packed
! triangular format
      integer j ! Column index
      integer j1 ! First column of block
      integer(long)  k ! Position in packed triangular format
      integer(long)  l ! Position in block packed triangular format
      integer lj ! Block length
      lj = n
      k = la + 1 - (n*(n+1_long))/2
      l = -lj
      do j1 =  1, np, nb
         do j = j1, min(j1+nb-1,np)
! Move column
            l = l + lj + 1
            a(l-j+j1:l-1) = 0.0_wp
            call dcopy(n+1-j,a(k),1,a(l),1)
            k = k + n+1-j
         end do
         lj = lj - nb
      end do
      kkn = l
      kdiag = k
      if (np<n) then
! To avoid traps on undefined values in calls of gemm from update2, set zeros
! in the part of a that is used as the first buffer and corresponds to
! unwanted values that are nevertheless passed to gemm.
        k = kkn + n+1-np
        do j = 1,min(n-p,nb/2)
          a(k:k+j-2) = 0.0_wp
          k = k + n-p
        end do
      end if
   end subroutine to_block

   subroutine from_block
! Rearrange columns 1 to q of a block packed triangular matrix so that each
! block column consists of the diagonal block packed by columns followed by
! the off-diagonal part as a rectangular matrix held by columns, and rearrange
! columns q+1 to np to packed triangular format.
      integer j ! Column index
      integer jb ! Number of columns in the block
      integer j1 ! First column of block
      integer(long) k ! Position in packed triangular format
      integer(long) k1 ! Position of the start of the block column
      integer(long) l ! Position in block packed triangular format
      integer(long) l1 ! New position in block packed triangular format
      integer lj ! Block length

      lj = n
      k1 =  1
      l1 =  1
      do j1 =  1, q, nb
         jb = min(nb,q-j1+1)
! Copy diagonal block to buf
         k =  1
         l = k1
         do j = 1, jb
            call dcopy(jb+1-j,a(l),1,buf(k),1)
            k = k + jb+1-j
            l = l + lj + 1
         end do
! Copy off-diagonal block to buf
         l = k1 + jb
         do j = 1, jb
            call dcopy(lj-jb,a(l),1,buf(k),1)
            k = k + lj-jb
            l = l + lj
         end do
         k1 = l - jb
! Copy the whole block column back
         k = 1
         do j = 1, jb
            call dcopy(lj+1-j,buf(k),1,a(l1),1)
            k = k+lj+1-j
            l1 = l1+lj+1-j
         end do
         lj = lj - nb
      end do
      ll = l1 - 1

! Rearrange columns q+1 to np to packed triangular format
      m = nb*(np/nb)
      q1 = 1 + nb*((q-1_long)/nb)
      lj = n - m
      l = kkn
      k = kdiag
      do j1 = m+1, q1, -nb
         do j = min(j1+nb-1,np), max(j1,q+1), -1
! Move column
            k = k - (n+1-j)
            call dcopy(n+1-j,a(l),1,a(k),1)
            l = l - lj - 1
         end do
         lj = lj + nb
       end do
   end subroutine from_block

   subroutine swap_cols(j1,j2)
   integer j1,j2
! Swap columns. It may be assumed that 1<=j1<j2<=p.

   integer d1,d2 ! Positions of the diagonal entries for columns j1,j2
   integer j ! Index of the leading column of the block
   integer jb ! Size of block
   integer(long) k1 ! Start of current part of column j1
   integer(long) k2 ! Start of current part of column j2
   integer l !
   integer lb ! No. rows in block column
   real(wp) temp

   j = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = j

! Swap columns of buffer
   call dswap(q-q1+1,buf(j1),n,buf(j2),n)

   k1 = j1
   k2 = j2
   lb = n

! Swap rows
   do j = 1, j1, nb
      jb = min(nb,n-j+1)
      l = min(nb,j1-j)
      if (l>0) call dswap(l,a(k1),lb,a(k2),lb)
      k1 = k1 + l*int(lb,long) - jb
      k2 = k2 + l*int(lb,long) - jb
      lb = lb - nb
   end do

   j = j - nb
   lb = lb + nb
   d1 = k1 + jb
   k1 = d1 + 1
   k2 = k2 + lb + jb

! Swap columns with rows
   jb = min(nb,j2-j1-1)
   do j = j, j2, nb
      l = min(jb,j+nb-1-j1,j2-j)
      if (l>0) call dswap(l,a(k1),1,a(k2),lb)
      k1 = k1 + l
      k2 = k2 + l*int(lb,long) - nb
      lb = lb - nb
   end do

   j = j - nb
   lb = lb + nb
   d2 = k2 + nb
   k2 = d2 + 1
   k1 = k1 + 1

! Swap the diagonals
   temp = a(d1)
   a(d1) = a(d2)
   a(d2) = temp

! Swap columns
   if (n>j2) call dswap(n-j2,a(k1),1,a(k2),1)

   end subroutine swap_cols
end subroutine ma64_factor_double



  subroutine ma64_solveL1_double(n,q,nb,b,flag,a,ll)
! Solves                ( L   ) x = b
!                       ( M  I)
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer(long) k ! Position in array
    integer lj ! Length of block

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return
    k = 1
    lj = n
    do j = 1, q, nb
      jb = min(nb,q-j+1)
      call dtpsv('L','N','U',jb,a(k),b(j),1)
      k = k + (jb*(jb+1_long))/2
      lj = lj - jb
      if (lj>0) call dgemv('N',lj,jb,-unity,a(k),lj,b(j),1,unity,b(j+jb),1)
      k = k + jb*lj
    end do

  end subroutine ma64_solveL1_double

  subroutine ma64_solveD1_double(n,q,b,flag,d)
! Solves                ( D   ) x = b
!                       (   I )
! where D is block diagonal of size q with blocks of size 1 and 2.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

! .. Locals ..
    integer j ! Column index
    real (wp) temp

    flag = 0
    if (n < 0) then
       flag = -1
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

! Apply operations associated with D
    j = 1
    do
      if (j>q) exit
      if (d(2*j)==0.0_wp) then
! 1x1 pivot
        b(j) = b(j)*d(2*j-1)
        j = j + 1
      else
! 2x2 pivot
        temp = b(j)*d(2*j-1) + b(j+1)*d(2*j)
        b(j+1) = b(j)*d(2*j) + b(j+1)*d(2*j+1)
        b(j) = temp
        j = j + 2
      end if
    end do

  end subroutine ma64_solveD1_double

  subroutine ma64_solveLT1_double(n,q,nb,b,flag,a,ll)
! Solves                ( L' M') x = b
!                       (    I )
! where  L is unit lower triangular of size q.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer j1 ! First column of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer m !

    flag = 0
! Back substitution
    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return
    m = nb*((q-1)/nb)
    k = 1 + ll
    j1 = 1 + m
    lj = n - q
    do j = j1, 1, -nb
      jb = min(nb,q-j+1)
      k = k - lj*jb
      if (lj>0) call dgemv('T',lj,jb,-unity,a(k),lj,b(j+jb),1,unity,b(j),1)
      k = k - (jb*(jb+1_long))/2
      call dtpsv('L','T','U',jb,a(k),b(j),1)
      lj = lj + jb
    end do

  end subroutine ma64_solveLT1_double

  subroutine ma64_solveDLT1_double(n,q,nb,b,flag,a,ll,d)
! Solves       ( D   ) ( L' M') x = b
!              (   I ) (    I )
! where D is block diagonal of size q with blocks of size 1 and 2
! and L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    real (wp), intent (inout) :: b(n) ! b(n) holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

    if (nb <= 1) then
       flag = -4
       return
    end if
    call ma64_solveD1_double(n,q,b,flag,d)
    call ma64_solveLT1_double(n,q,nb,b,flag,a,ll)

  end subroutine ma64_solveDLT1_double

  subroutine ma64_solveL2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)
! Solves                ( L   ) x = b
!                       ( M  I)
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer r ! Index for rhs

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

     k = 1
     lj = n
     do j = 1, q, nb
        jb = min(nb,q-j+1)
        do r = 1,nrhs
           call dtpsv('L','N','U',jb,a(k),b(j,r),1)
        end do
        k = k + (jb*(jb+1_long))/2
        lj = lj - jb
        if (lj>0) call dgemm &
           ('N','N',lj,nrhs,jb,-unity,a(k),lj,b(j,1),ldb,unity,b(j+jb,1),ldb)
        k = k + jb*int(lj,long)
     end do

  end subroutine ma64_solveL2_double

  subroutine ma64_solveD2_double(n,q,nrhs,b,ldb,flag,d)
! Solves                ( D   ) x = b
!                       (   I )
! where D is block diagonal of size q with blocks of size 1 and 2.

   integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

! .. Locals ..
    integer j ! Column index

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q==0) return

    do j = 1, nrhs
       call ma64_solveD1_double(n,q,b(1,j),flag,d)
    end do

  end subroutine ma64_solveD2_double

  subroutine ma64_solveLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)
! Solves                ( L' M') x = b
!                       (    I )
! where  L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of b
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )

! .. Locals ..
    integer j ! Column index
    integer jb ! Width of block
    integer j1 ! First column of block
    integer(long) k ! Position in array
    integer lj ! Length of block
    integer m !
    integer r ! Index for rhs

    flag = 0
    if (n < 0) then
       flag = -1
    else if (nb <= 1) then
       flag = -4
    else if (nrhs < 0) then
       flag = -5
    else if (ldb < n) then
       flag = -6
    else if (q < 0) then
       flag = -8
    else if (q > n) then
       flag = -9
    end if
    if (flag/=0 .or. q<1) return

! Back substitution
     m = nb*((q-1)/nb)
     k = 1 + ll
     j1 = 1 + m
     lj = n - q
     do j = j1, 1, -nb
        jb = min(nb,q-j+1)
        k = k - jb*int(lj,long)
        if (lj>0) call dgemm &
            ('T','N',jb,nrhs,lj,-unity,a(k),lj,b(j+jb,1),ldb,unity,b(j,1),ldb)
        k = k - (jb*(jb+1_long))/2
        do r = 1, nrhs
           call dtpsv('L','T','U',jb,a(k),b(j,r),1)
        end do
        lj = lj + jb
     end do

  end subroutine ma64_solveLT2_double

  subroutine ma64_solveDLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll,d)
! Solves       ( D   ) ( L' M') X = B
!              (   I ) (    I )
! where D is block diagonal of size q with blocks of size 1 and 2
! and L is unit lower triangular of size q.

    integer, intent (in) :: n ! Matrix order
    integer, intent (in) :: q ! Order of L (and D)
    integer, intent (in) :: nb ! The block size used for the blocked hybrid
!                                format.
    integer, intent (in) :: nrhs ! Number of right-hand sides
    integer, intent (in) :: ldb ! Leading extent of B
    real (wp), intent (inout) :: b(ldb,nrhs) ! Holds the right-hand
!                           sides on entry and is overwritten by the solution.
    integer, intent (out) :: flag
    integer(long), intent (in) :: ll ! Size of a
    real (wp), intent (in) :: a( ll )
!                               Holds ( L ) packed by columns.
!                                     ( M )
    real (wp), intent (in) :: d(2*q) !  Holds the inverse of D,
!              except that zero diagonal blocks (pivots) are not inverted.
!              Diagonal entries are in d(1:2*q-1:2) and entries below
!              the diagonal are in d(2:2*q-2:2).

    if (nb <= 1) then
       flag = -4
       return
    end if
    call ma64_solveD2_double(n,q,nrhs,b,ldb,flag,d)
    if (flag/=0) return
    call ma64_solveLT2_double(n,q,nb,nrhs,b,ldb,flag,a,ll)

  end subroutine ma64_solveDLT2_double


end module hsl_ma64_double


! COPYRIGHT (c) 2007 Council for the Central Laboratory
!               of the Research Councils
! Original date 5 October 2007. Version 1.0.0.
! 13 November 2007. Version 1.1.0 Subroutines of01*copy moved into the module.
!   and removed the separate data module.
! 22 November 2007. Version 2.0.0 Logical argument active replaced by
!    integer(long) argument inactive.
! 7 December 2007. Version 3.0.0 Argument inactive removed from OF01_read
!    and argument retain replaced by discard.
! 29 July 2009. Version 3.1.0 Minor changes made to avoid references to 
!    undefined variables.
! 17 December 2009. Version 3.2.0 Minor changes stay within 80 characters
!    per line. 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! To generate another version 
! (1) Globally change  _double
! (2) Globally change dcopy 
! (3) Change wp = kind(0.0d0) 
! (4) Globally change the type real(wp)

module hsl_of01_double
  implicit none
  private

  public of01_data, of01_data_private
  public of01_initialize, of01_open, of01_close, of01_read, of01_write, &
    of01_end


  !!! Parameters !!!
  integer, parameter  :: default_lpage = 2**12 ! Default page length.

  integer, parameter  :: default_maxfiles = 10 ! Maximum number of open files.

  integer, parameter  :: default_npage = 1600 ! Default number of pages
!          in the buffer.

  integer, parameter  :: long = selected_int_kind(18) ! Long integer.

  integer, parameter  :: default_file_size = 2**21 ! Default
!          target length of each file.

  integer, parameter  :: nsup = 2 ! Initial size of the array
!          data%private%filename.

  integer, parameter  :: wp = kind(0.0d0) ! Defines data type (single/double).

  integer, parameter  :: ihash = 3 ! Constant used in the hashing function.

  integer, parameter  :: maxpath = 400 ! Max. length of path name.

  integer, parameter  :: maxname = 400 ! Max. length of superfile name.

  integer, parameter  :: maxlen = maxpath+maxname+10 ! Max. length of file name
!          There has to be a limit on the lengths of path and file names because
!          Fortran 95 requires the file name in an OPEN statement to be a scalar
!          character variable.

  real(wp), parameter :: zero = 0.0_wp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type of01_data_private

! Note: Age is measured since last reference.
    real(wp), allocatable :: buffer(:,:) ! of shape (lpage,npage)

    logical,  allocatable :: differs(:) ! of size npage.
!             differs(I)=.true. if the file version of page I
!             of the buffer is different from the buffer version.

    character(maxname), allocatable :: filename(:)
!             Set to hold a copy of the superfile names given to OF01_open.

    integer(long), allocatable :: first(:) ! of size npage. first(K) is the
!            first page with hash code K or zero if there are none.

    integer :: free !  Start of linked list of free file indices.

    integer(long), allocatable :: left(:),right(:) ! of size maxfiles.
!            A sequence of discarded entries on a superfile is recorded as
!            left(superfile):right(superfile).

    integer(long), allocatable :: highest(:) ! of size maxfiles.
!            highest(I) holds:
!            -1 after call of of01_initialize,
!            -2 after call of of01_close for file I, or
!            >0 highest page number read or written on file I.

    integer, allocatable :: index(:) ! of size npage. index(I) contains the
!            index of the (super)file associated with page I of the buffer.

    integer :: iolength ! iolength of a page

    integer :: maxfiles = 0 ! Maximum number of open files.

    integer, allocatable :: name(:) ! of size maxfiles. For a primary file,
!            name(I) holds the position of the superfile name in filename.

    integer(long), allocatable :: next(:) ! of size npage.
!            next(I) is the next page with same hash code
!            or zero for last in list. Otherwise undefined.

    integer :: nfiles !  Number of file indices in use.

    integer(long) :: nrec ! Number of records in each file.

    integer, allocatable :: nextfile(:) ! of size maxfiles. Used for the
!            linked list of free file indices and linked lists of indices of
!            files in a superfile.  nextfile(I) is the next file index in
!            its list or zero if there are none.

    integer(long), allocatable :: older(:) ! of size npage. older(I) is the
!            next older page to page I, or the youngest page if I is oldest.

    integer(long), allocatable :: page(:) ! of size npage. page(I) is the
!            file address (page number) of page I in the buffer.

    integer(long), allocatable :: page_list(:) ! of size npage.
!            list of pages for this read or write
!            that were in the buffer at the time of call.

    character(maxpath), allocatable :: path(:)
!            The paths given to OF01_initialize

    integer(long), allocatable :: prev(:) ! of size npage. prev(I) is the
!            previous page with same hash code or -(hash code)
!            for first in list. Otherwise undefined.

    integer, allocatable :: unit(:) ! of size maxfiles. unit(I) holds the
!            unit for file I or zero if not in use.

    integer(long), allocatable :: younger(:) ! of size npage.
!            younger(I) is the next younger page to page I, or
!            the oldest page if I is the youngest.

    integer :: youngest ! most recently referenced page in the buffer

  end type of01_data_private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type of01_data

    integer :: entry = 0 ! Index of the entry to the package.

    integer :: iostat    ! Fortran iostat parameter.
    integer :: lpage     ! length of each page of the buffer, that is, the
!            number of scalar variables in each page. The default is 1024.

    integer(long) :: ncall_read  ! number of calls to OF01_read.

    integer(long) :: ncall_write ! number of calls to OF01_write.

    integer(long) :: nio_read ! number of records read by
!            OF01_read and OF01_write.

    integer(long) :: nio_write ! number of records written by
!            OF01_read and OF01_write.

    integer(long) :: npage ! number of pages in the in-core buffer.
!            The default value is 20. It is a long integer, but
!            we do not anticipate it being very large.

    integer(long) :: file_size ! target length of each file.
!             The default value is default_file_size.

    integer(long) :: nwd_read ! number of scalars read by OF01_read.

    integer(long) :: nwd_write ! number of scalars written by OF01_write.

    type (of01_data_private) :: private
! In Fortran 2003, this should be
!   type (of01_data_private), private :: private

    integer :: stat ! Fortran stat parameter.

  end type of01_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface of01_initialize
    module procedure of01_initialize_double
  end interface
  interface of01_open
    module procedure of01_open_double
  end interface
  interface of01_close
    module procedure of01_close_double
  end interface
  interface of01_read
    module procedure of01_read_double
  end interface
  interface of01_write
    module procedure of01_write_double
  end interface
  interface of01_end
    module procedure of01_end_double
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_initialize_double &
                            (iflag,data,path,file_size,lpage,npage,lp)

! This subroutine must be called once to initialize
! a structure of derived type of01_data
! and to allocate and initialize its array components.

    integer, intent (out) :: iflag ! Value 0 on successful return.
!            The only possible negative values are:
!            -1 Allocation error. The STAT parameter returned in data%stat.
!            -8 Deallocation error. The STAT parameter returned in data%stat.
!           -16 Character length of path too great.

    type (of01_data) , intent (inout) :: data

    character(*), optional, intent (in) :: path(:) ! Paths for the files.
!            If absent, (\ '' \) is used.

    integer(long), optional, intent (in) :: file_size ! If present and positive,
!            the target length of a file. Default value
!            used if absent or present and not positive.

    integer, optional, intent (in) :: lpage ! If present and positive, the
!            size of each page in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: npage ! If present and positive, the
!            number of pages in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: lp ! If present and not negative, the
!            unit number for diagnostic messages. 6 is
!            used if absent or present and negative.

! local variables
    integer :: i !  do loop variable
    integer :: size_path ! size of array path

!!!!!!!!!!!!!!!!!!!!!!!
! Initialise
    iflag = 0
    data%entry       = 1
    data%ncall_read  = 0
    data%ncall_write = 0
    data%nio_read    = 0
    data%nio_write   = 0
    data%nwd_read    = 0
    data%nwd_write   = 0

! If napge/lpage supplied, use instead of the default values
    data%npage = default_npage
    if (present(npage)) then
      if (npage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%npage = npage
    end if

    data%lpage = default_lpage
    if (present(lpage)) then
      if (lpage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%lpage = lpage
    end if

! Find record length
    if (allocated(data%private%buffer)) then
       deallocate (data%private%buffer,stat=data%stat)
    end if
    allocate (data%private%buffer(data%lpage,1),stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -1; call print_iflag(data,iflag,lp); return
    end if
    data%private%buffer(:,1) = 0
    inquire (iolength=data%private%iolength) data%private%buffer(:,1)
    deallocate (data%private%buffer,stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -8; call print_iflag(data,iflag,lp); return
    end if

    size_path = 1
    if (present(path)) then
       size_path = size(path)
       if (len(path)>maxpath) then
         iflag = -16; call print_iflag(data,iflag,lp); return
       end if
    end if

    data%file_size = default_file_size
    if (present(file_size)) then
      if (file_size < data%lpage) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%file_size = file_size
    end if
    data%private%nrec = data%file_size/data%lpage
    data%file_size = data%private%nrec*data%lpage

! Allocate private arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,   data%private%nextfile, &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%path,      data%private%name,     &
                  data%private%left,      data%private%right,    &
                  data%private%filename, stat=data%stat)
       if (data%stat /= 0) then
         iflag = -8; call print_iflag(data,iflag,lp); return
       end if
    end if
    data%private%maxfiles = default_maxfiles
    allocate (data%private%highest(data%private%maxfiles),  &
              data%private%nextfile(data%private%maxfiles), &
              data%private%index(1:data%npage),             &
              data%private%older(1:data%npage),             &
              data%private%younger(1:data%npage),           &
              data%private%differs(1:data%npage),           &
              data%private%prev(1:data%npage),              &
              data%private%next(1:data%npage),              &
              data%private%first(1:data%npage),             &
              data%private%page(1:data%npage),              &
              data%private%page_list(1:data%npage),         &
              data%private%unit(data%private%maxfiles),     &
              data%private%path(size_path),                 &
              data%private%name(data%private%maxfiles),     &
              data%private%left(data%private%maxfiles),     &
              data%private%right(data%private%maxfiles),    &
              data%private%filename(nsup),                  &
              data%private%buffer(data%lpage,data%npage), stat=data%stat)
    if (data%stat /= 0) then
      data%private%maxfiles = 0
      iflag = -1; call print_iflag(data,iflag,lp); return
    end if

! Initialise arrays
    data%private%free = 0
    data%private%nfiles = 0
    data%private%filename(:) = ''
    do i = 1, data%npage
      data%private%index(i) = -1
      data%private%older(i) = i + 1
      data%private%younger(i) = i - 1
      data%private%page(i) = 0
      data%private%first(i) = 0
      data%private%next(i) = 0
      data%private%prev(i) = 0
      data%private%differs(i) = .false.
    end do
    data%private%youngest = 1
    data%private%younger(1) = data%npage
    data%private%older(data%npage) = 1
    if (present(path)) then
       data%private%path(:) = path(:)
    else
       data%private%path(:) = (/ '' /)
    end if
    data%private%name(:) = 0
    data%private%buffer(:,1) = 0

  end subroutine of01_initialize_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_open_double(filename,ifile,iflag,data,lenw,lp)

! This subroutine must be called for each superfile that is to be
! accessed through HSL_OF01. It gives the superfile an index and opens
! its files.

    character(len=*), intent (in) :: filename ! Name of the superfile.

    integer, intent (out) :: ifile ! Index of the superfile.

    integer, intent (out) ::  iflag ! Value 0 on successful return.
!             A negative value is associated with an error message
!               on unit lp. Possible negative values are:
!            -1  Allocation error. The STAT parameter returned in data%stat.
!            -5  Error in inquire statement.
!            -7  Error in open statement.
!            -8  Deallocation error. The STAT parameterreturned in data%stat.
!            -11 lenw > 0, but not enough files exist.
!            -12 File filename already exists, but lenw is not present
!                or lenw <=  0.
!            -13 filename is too long.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer(long), optional, intent (in) :: lenw ! length in pages of the
!            part of the file that has been written and is not
!            regarded as having been overwritten by zeros.
!            Pages beyond this are regarded as containing zeros.
!            If not present or lenw <=  0, a new file is opened.

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.
!            If lp < 0, messages are suppressed.

! Local variables
    character(10) :: ci ! Filename extension
    integer :: i ! do loop variable
    integer :: k ! number of secondary files
    integer :: l ! file index
    integer(long) :: lenw_copy ! no. of pages in the virtual array
    integer :: nout ! unit for error messages
    integer :: m ! Previous file index in linked list
    character(maxname), allocatable :: temp(:)

    iflag = 0
    nout = 6
    if (present(lp)) nout = lp
    data%entry = 2

    if (len(filename) > maxname) then
       iflag = -13; go to 100
    end if

! Find the length of the virtual array
    lenw_copy = 0
    if (present(lenw)) lenw_copy = max(0_long,lenw)

! Find number of secondary files
    k = max(0_long,(lenw_copy-1)/data%private%nrec)

! Find suitable indices and units for all the files and open them
    m = 0
    ci=''
    do i = 0,k
      if (lenw_copy > 0) then
        call open_old (data,l,trim(filename)//adjustl(ci),iflag)
      else
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
      end if
      if (iflag /= 0) go to 100
      if (i==0) ifile = l
      data%private%highest(l) = max(0_long,lenw_copy)
      data%private%left(l) = 1
      data%private%right(l) = 0
      if (m > 0) data%private%nextfile(m) = l
      m = l
      data%private%highest(l) = lenw_copy
      write(ci,'(i5)') i + 1
    end do

! Look for a place and store filename there
    k = size(data%private%filename)
    do i = 1, k
      if (data%private%filename(i)== '') then
          data%private%filename(i) = filename
          data%private%name(ifile) = i
          return
      end if
    end do
! Increase size of data%private%filename
    allocate(temp(k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    temp(:) = data%private%filename(:)
    deallocate(data%private%filename,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
    allocate(data%private%filename(2*k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    data%private%filename(1:k) = temp(:)
    data%private%filename(k+1:) = ''
    deallocate(temp,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
! Store filename
    data%private%filename(k+1) = filename
    data%private%name(ifile) = k+1
    return

100 call print_iflag(data,iflag,lp)

  end subroutine of01_open_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_read_double &
       (ifile,loc,n,read_array,iflag,data,lp,map,discard)

! This subroutine performs the read from a superfile.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be read.

    real(wp), intent (inout) :: read_array(*)
!            Array into which data from the file is read.
!            If map is present, file data is added into it thus:
!               read_array(map) = read_array(map) + file_data(1:n)

    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not open.
!            -5  Error in Fortran INQUIRE. The IOSTAT parameter is returned in
!                data%iostat.
!            -6  Error in Fortran READ. The IOSTAT parameter is returned in
!                data%iostat.
!            -7  Error in Fortran OPEN statement. The IOSTAT parameter
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

    integer, optional, intent (in) :: map(n) ! map array.

    logical, optional, intent (in) :: discard ! Whether data is to be discarded.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos  ! position within first page of first value
!           required.
    integer(long) :: left,right ! Range of discarded entries
    integer(long) :: ip  ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer
!           pages that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer       :: jfile ! file associated with a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!           required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer       :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: my_discard ! Whether the data is to be discarded.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    if (n <= 0) return
    data%entry = 3

! Check for errors in incoming data
    if (ifile < 1) then
       iflag = -9; go to 200
    else if (ifile > data%private%nfiles) then
       iflag = -4; go to 200
    else if (data%private%highest(ifile) < 0) then
       iflag = -4; go to 200
    else if (loc <= 0) then
       iflag = -3; go to 200
    end if

    my_discard = .false.
    if (present(discard)) then
      if (discard) then
        my_discard = .true.
        left = loc
        right = loc+n-1
        if (loc == data%private%right(ifile)+1) then
          left = data%private%left(ifile)
        else if (loc+n == data%private%left(ifile)) then
          right = data%private%right(ifile)
        end if
        data%private%left(ifile) = left
        data%private%right(ifile) = right
      end if
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_read = data%ncall_read + 1
    data%nwd_read = data%nwd_read + n

    lenl = 0
! Look for required pages that are in the buffer
    look1: do ip = first_page, last_page

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp) == ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
10    call read_from_buffer
      if (my_discard) then
        if ((ip-1)*lpage+1 >= left .and. ip*lpage+1 <= right ) &
          data%private%differs(jp) = .false.
      end if
   end do look1

! look for pages not yet treated and treat each of them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw.le.lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if
! the oldest page is overwritten and becomes the youngest
      jp = data%private%younger(data%private%youngest)
      data%private%youngest = jp

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile>=0) then
! write page from buffer to record page(jp) in the file
! if it has been altered (ie if differs(,jp) is set to 1).
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 200

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip > data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 200
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m > 0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .false.

! perform actual transfer to read_array from buf
      call read_from_buffer

    end do look2
!    call buffer
    return

200 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_from_buffer

      integer :: l1,l2 ! Temporaries
! Read page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip == last_page) l2 = last_pos
      if (.not. present(map)) then
        call dcopy(l2-l1+1,data%private%buffer(l1,jp),1,read_array(m),1)
      else
        call of01_mcopy(l2-l1+1,data%private%buffer(l1,jp),read_array,map(m))
      end if
    end subroutine read_from_buffer

  end subroutine of01_read_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_write_double &
             (ifile,loc,n,write_array,iflag,data,lp,inactive)

! This subroutine performs the writing to a superfile.

    integer, intent (in) :: ifile ! Index of the (super)file.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be written.

    real(wp), intent (in) :: write_array(*) ! array from which data
!            is written to the file.
    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not OF01_open.
!            -5  Error in INQUIRE.  The IOSTAT value
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.

    integer(long), optional, intent (in) :: inactive ! If present, it is
!            assumed that the data transferred are unlikely to be needed before
!            other data in the buffer. If inactive<loc, it is assumed that this
!            is true also for the data in the interval inactive:loc. If
!            inactive loc+n, it is assumed that this  is true also for the data
!            in the interval loc+n:inactive.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos !  position within first page of first value
!            required.
    integer(long) :: ip ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer pages
!            that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer :: jfile ! file index of a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!            required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: inactive_page ! Whether this page is fully inactive.
    logical       :: inactive_first ! Whether first page is fully inactive.
    logical       :: inactive_middle ! Whether middle pages are fully inactive.
    logical       :: inactive_last ! Whether last page is fully inactive.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    if (n <= 0) return
    data%entry = 4

! Check for errors in incoming data
    if (ifile < 1) then
      iflag = -9; go to 300
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 300
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 300
    else if (loc <= 0) then
      iflag = -3; go to 300
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

    if (present(inactive)) then
       if (first_page==last_page) then
          inactive_first = .false.
       else
          inactive_first = inactive<=1+(first_page-1)*lpage
          inactive_middle = .true.
          inactive_last = inactive>=last_page*lpage
       end if
    else
       inactive_middle = .false.
       inactive_first = .false.
       inactive_last = .false.
    end if

! If there is any overlap with the range of discarded entries,
! reset the range to be null
    if (loc>data%private%right(ifile)) then
    else if (loc+n-1<data%private%left(ifile)) then
    else
       data%private%left(ifile) = 1
       data%private%right(ifile) = 0
    end if

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_write = data%ncall_write + 1
    data%nwd_write = data%nwd_write + n

    lenl = 0
! Look for required pages in the buffer
    look1: do ip = first_page, last_page
      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if
!     if (inactive_page) write(10,'(i3,i9,3i5)') ifile,loc,n,ip,data%lpage

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp)==ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          call write_to_buffer
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
      call write_to_buffer
! If the page is fully inactive, make it the oldest
 10   if (inactive_page) then
         data%private%youngest = data%private%older(jp)
      end if
    end do look1

! look for pages not yet treated and treat them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw <= lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if

      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if

! the oldest page is overwritten
      jp = data%private%younger(data%private%youngest)
! if transfer is active or the page is partially active, page becomes
! the youngest
     if (.not.inactive_page) then
        data%private%youngest = jp
     end if

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile >= 0) then
! write page from buffer to record page(jp) in the file
! if it has been altered.
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 300

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip>data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 300
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m>0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .true.

! perform actual transfer to write_array from buf
      call write_to_buffer

    end do look2
!    call buffer
    return

300 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_to_buffer

      integer :: l1,l2 ! Temporaries

! Write page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip==last_page) l2 = last_pos
      call dcopy(l2-l1+1,write_array(m),1,data%private%buffer(l1,jp),1)
      data%private%differs(jp) = .true.

    end subroutine write_to_buffer
  end subroutine of01_write_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_close_double(ifile,lenw,num_file,iflag,data,lp,lkeep)

! This subroutine should be called when a superfile is finished with.
! Its data in the buffer are unloaded.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (out) :: lenw ! length in pages of the part of
!            the file that has been written.
    integer, intent (out) :: num_file ! total number of additional files
!            that have been opened. They are named
!            filenamej, j = 1, ..., num_file.
    integer, intent (out) :: iflag ! iflag = 0 for a successful return.
!            A negative value is associated with an error
!            message on unit lp. Possible negative values:
!            -4  Attempt to access a file that is not OF01_open.
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -9  ifile out of range.
!            -14 Error in CLOSE statement. The IOSTAT value
!                is returned in data%iostat.
!            -15 Error in WRITE.  The IOSTAT value
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0 , messages are suppressed.

    logical, optional, intent (in) :: lkeep ! If present and set to .FALSE.,
!            the files are deleted on being closed.
!            Otherwise, the files are kept.

! Local variables
    integer(long) :: i !  do loop variable
    integer :: j       !  do loop variable
    integer :: jfile   !  index of a file
    integer(long) :: l ! page number
    integer :: lpage   ! length of each page (= data%lpage)
    integer(long) :: m ! page number
    integer(long) :: oldest ! oldest page
    integer(long) :: oldest2 ! second oldest page
    character (maxlen) name  ! File name
    character(len=6) :: status ! File status on closing

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 5
! Check for errors
    if (ifile < 1) then
      iflag = -9; go to 400
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 400
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 400
    end if

    status = 'keep'
    if (present(lkeep)) then
      if (.not. lkeep) status = 'delete'
    end if

! search the pages of the buffer dealing with those associated
! with ifile.
    search:do i = 1, data%npage
! data%private%index(i) is index of the (super)file associated
! with page i of buffer
      if (ifile /= data%private%index(i)) cycle
      data%private%index(i) = -1
      if (data%private%differs(i)) then
! File version of page i is not the same as the buffer version so
! write page i to file
        data%private%differs(i) = .false.
        lpage = data%lpage
        if (status == 'keep') then
          call page_write(data, &
          data%private%buffer(:,i),ifile,data%private%page(i),iflag)
          if (iflag /= 0) then
            call print_iflag(data,iflag,lp)
            return
          end if
        end if
      end if

! remove page from hash list
      l = data%private%prev(i) ! Previous
      m = data%private%next(i) ! Next
      if (m>0) data%private%prev(m) = l
      if (l>0) then
        data%private%next(l) = m
      else if (l < 0) then
        l = -l  ! i is first in the list with code l      
        data%private%first(l) = m
      end if

! youngest is the most recently referenced page in the buffer.
! younger(youngest) is the oldest page in the buffer.
      oldest = data%private%younger(data%private%youngest)
! if page i is not already the oldest, make it so.
      if (oldest==i) cycle
! remove page from its present position
      l = data%private%younger(i)
      m = data%private%older(i)
! l is next youngest to page I and m is the next oldest
      data%private%older(l) = m
      data%private%younger(m) = l

! insert page i as the oldest page
! older(oldest) is the youngest page
      data%private%youngest = data%private%older(oldest)

! Page i becomes the oldest page. The page that was the oldest
! becomes the second (ie next) oldest page
      oldest2 = oldest
      oldest = i
      data%private%older(oldest2) = oldest
      data%private%younger(oldest) = oldest2
      data%private%older(oldest) = data%private%youngest
      data%private%younger(data%private%youngest) = oldest
    end do search

    lenw = data%private%highest(ifile)
    if (lenw <= 0 ) status = 'delete'
    num_file = (lenw-1)/data%private%nrec

! remove file from list of declared files.
! close files
    inquire(unit=data%private%unit(ifile),name=name,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -5; go to 400
    end if
    close (unit=data%private%unit(ifile),status=status,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -14; go to 400
    end if
    data%private%unit(ifile) = 0
    data%private%highest(ifile) = -1
    jfile = ifile
    if (num_file > 0) then
      do j = 1, num_file
        m = data%private%nextfile(jfile)
        inquire(unit=data%private%unit(m),name=name,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -5; go to 400
        end if
        close (unit=data%private%unit(m),status=status,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -14; go to 400
        end if
        data%private%unit(m) = 0
        data%private%highest(m) = -1
        jfile = m
      end do
    end if
! Put the freed indices at the front of the list list of free indices
    data%private%nextfile(jfile) = data%private%free
    data%private%free = ifile

! Clear the name for reuse
    data%private%filename(data%private%name(ifile))=''

    return

400 call print_iflag(data,iflag,lp)

  end subroutine of01_close_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_end_double(iflag,data,lp)

! This subroutine deallocates the private array components of data.

    integer, intent (out) :: iflag
!             0  successful return.
!            -8  Error in Fortran DEALLOCATE statement. The STAT
!                value is returned in data%stat.
!            -10 of01_close has not been called for a superfile.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0, messages are suppressed.

! Local variables
    integer :: i

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 6

    do i = 1, data%private%nfiles
      if (data%private%highest(i) >= 0) then
        iflag = -10; call print_iflag(data,iflag,lp)
        return
      end if
    end do

! Deallocate arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,  data%private%nextfile,  &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%buffer,    data%private%path,     &
                  data%private%left,      data%private%right,    &
                  data%private%name,      data%private%filename, &
                  stat=data%stat)
      data%private%maxfiles = 0
      if (data%stat /= 0) then
        iflag = -8; call print_iflag(data,iflag,lp)
      end if
    end if

  end subroutine of01_end_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_write(data,buffer,ifile,locp,iflag)

! This subroutine transfers a page of the buffer to a file. If the
! page is within the existing virtual array then a simple write
! is performed. To write a page beyond the present end may require
! writing pages full of zeros in front of its position.

    type (of01_data), intent (inout) :: data

    real(wp), intent (inout) :: buffer(data%lpage) ! array from which
!             the transfer is made.
    integer, intent (in) :: ifile ! the index of the primary file.

    integer(long), intent (in) :: locp ! page number in the file.

    integer, intent (out) :: iflag ! error flag. Possible negative values:
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -7  Error in OPEN. IOSTAT valueis returned in data%iostat.
!            -15 Error in WRITE. IOSTAT value is returned in data%iostat.

! Local variables
    integer(long) :: hpage ! highest page read or written for ifile
    integer(long) :: l ! do loop variable

    hpage = data%private%highest(ifile)
    if (locp <= hpage+1) then
! Write straight to file
      call actual_page_write(locp)
      if (iflag /= 0) return
    else
! Fill the gap with zero page
! first write out buffer to page hpage+1
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
      buffer(:) = zero
! write pages of zeros
      do l = hpage+2, locp-1
        call actual_page_write(l)
        if (iflag /= 0) return
      end do
! read buffer back in from page hpage+1
      call page_read(data,buffer,ifile,hpage+1,iflag)
      if (iflag /= 0) return
! Now write buffer to the required location
      call actual_page_write(locp)
      if (iflag /= 0) return
! Finally write zeros to hpage+1
      buffer(:) = zero
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
    end if

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine actual_page_write(locp)

! Perform an actual i/o write of a page to a superfile. If necessary, open
! intermediate secondary files.

    integer(long), intent (in) :: locp ! page index in the superfile.

! Local variables
    character(10) :: ci ! Filename extension.
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or >0 for k-th secondary file.
    integer l ! index of the file that is written.
    integer m ! Previous file in linked list.
    character(len=maxlen) :: filename ! name of the primary file.

    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec
    l = ifile
    do i = 1, k
      m = l
      l = data%private%nextfile(l)
      if (l == 0) then
! open new file (filename1, filename2, etc)
        filename = data%private%filename(data%private%name(ifile))
        write(ci,'(i5)')i
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
        if (iflag /= 0)  return

! Add to linked list for ifile
        if (m > 0) data%private%nextfile(m) = l
        m = l
      end if
    end do

    data%nio_write = data%nio_write + 1
    write (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
!    write(*,*) 'write record',new_locp,data%iostat
    if (data%iostat /= 0) then
      iflag = -15; return
    end if
    data%private%highest(ifile) = max(data%private%highest(ifile),locp)

  end subroutine actual_page_write

  end subroutine page_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_read(data,buffer,ifile,locp,iflag)

! This subroutine performs an actual i/o read of a page from a superfile

    type (of01_data), intent (inout) :: data

    real(wp),intent(out) :: buffer(data%lpage) ! array to which page is read

    integer, intent (in) :: ifile ! index of the superfile.

    integer(long), intent (in) :: locp ! page index in the superfile.

    integer, intent (out) :: iflag ! Error flag. Possible negative value:
!            -6 error in read statement

! Local variables
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or > 0 for k-th secondary file.
    integer l ! index of the file that is read.

!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec

! Find index of the kth linked file
    l = ifile
    do i = 1, k
      l = data%private%nextfile(l)
    end do

    data%nio_read = data%nio_read + 1
    read (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
    if (data%iostat /= 0) iflag = -6

  end subroutine page_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_new (data,ifile,filename,iflag)

! Open a new file

    type (of01_data) :: data

    integer, intent(out) :: ifile ! File index.

    character (*) :: filename!    ! File name.

    integer, intent (inout) :: iflag ! Error flag. Possible nonzero value:
!            -7 Error in OPEN.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: j ! Do index.
    integer :: unit ! Unit.
    logical :: open ! open variable for inquire.

! Check file does not already exist
    do i = 1, size(data%private%path)
       inquire (file=trim(data%private%path(i))//filename,exist=ex,&
          iostat=data%iostat)
       if (data%iostat /= 0) then
          iflag = -5; return
       else if (ex) then
          iflag = -12; data%iostat=i; return
       end if
    end do

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

! Open the file on the first available path
    outer: do i = 1, size(data%private%path)
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       if (data%iostat /= 0) then
         close(unit,iostat=data%iostat,status='delete')
         cycle
       end if

       if (size(data%private%path) == 1) return
! Check that there is room for the file
       do j = 1, data%private%nrec
         write(unit,rec=j,iostat=data%iostat)data%private%buffer(:,1)
         if (data%iostat /= 0) then
           close(unit,iostat=data%iostat,status='delete')
           cycle outer
         end if
       end do
       close(unit,iostat=data%iostat)
       if (data%iostat /= 0) then
         inquire(unit,opened=open)
         if (.not.open) &
         open (unit,file=trim(data%private%path(i))//filename, &
             access='direct',iostat=data%iostat,recl=data%private%iolength)
         close(unit,iostat=data%iostat,status='delete')
         cycle outer
       end if
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       return
    end do outer
    iflag = -17

  end subroutine open_new

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_old (data,ifile,filename,iflag)

! Open an old file.

    type (of01_data) :: data

    integer, intent(out) :: ifile ! File index.
    character (*) :: filename     ! File name.

    integer, intent (inout) :: iflag ! Possible nonzero values:
!            -5  Error in inquire statement
!            -7  Error in OPEN.  The IOSTAT value is returned in data%iostat.
!            -11 File does not exist

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: unit ! Unit.

! Find the file.
    do i = 1, size(data%private%path)
      inquire (file=trim(data%private%path(i))//filename,exist=ex,&
             iostat=data%iostat)
      if (data%stat /= 0) then
        iflag = -5; return
      end if
      if (ex) exit
    end do
    if (.not. ex) then
! File not found
      iflag = -11; return
    end if

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

    open (unit,file=trim(data%private%path(i))//filename, access='direct', &
          iostat=data%iostat, recl=data%private%iolength)
    if (data%iostat /= 0) then
      iflag = -7; return
    end if

  end subroutine open_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_unit(data,unit,ifile,iflag)

! Find a free unit and a free file index.

    type (of01_data) :: data
    integer, intent (out) :: unit ! Free unit.

    integer, intent (out) :: ifile ! Free file.

    integer, intent (inout) :: iflag ! Possible nonzero value:
!           -5 Error in INQUIRE.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: open
    integer :: oldsize

! Find free unit
    do unit = 12, huge(0)
      inquire (unit=unit,iostat=data%iostat,opened=open)
      if (data%iostat /= 0) then
! Error in inquire
        ifile = -1
        iflag = -5; return
      end if
      if (.not. open) exit
    end do

! Find free file index
    ifile = data%private%free
    if (ifile > 0) then
      data%private%free = data%private%nextfile(ifile)
    else
      ifile =  data%private%nfiles+1
      if (ifile>data%private%maxfiles) then
! Reallocate files
        data%private%maxfiles = 0
        call reallocate_long(data%private%highest)
        if (iflag < 0) return
        call reallocate_long(data%private%left)
        if (iflag < 0) return
        call reallocate_long(data%private%right)
        if (iflag < 0) return
        call reallocate(data%private%nextfile)
        if (iflag < 0) return
        oldsize = size(data%private%name)
        call reallocate(data%private%name)
        data%private%name(oldsize+1:) = 0
        if (iflag < 0) return
        call reallocate(data%private%unit)
        if (iflag < 0) return
        data%private%maxfiles = size(data%private%unit)
      end if
      data%private%nfiles =  ifile
    end if
    data%private%nextfile(ifile) = 0
    data%private%highest(ifile) = -1
    data%private%left(ifile) = 1
    data%private%right(ifile) = 0
    data%private%unit(ifile) = unit

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    subroutine reallocate(array)

      integer, allocatable :: array(:)
      integer, allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)
      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reallocate_long(array)

      integer(long), allocatable :: array(:)

      integer(long), allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)

      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate_long

  end subroutine find_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_iflag(data,iflag,lp)

    type (of01_data) :: data

    integer, intent (in) :: iflag ! Error flag.

    integer, optional :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

! Local variables
    integer :: i, nout
    character(10), parameter :: names(6) = (/'initialize','open      ',&
                    'read      ','write     ','close     ','end       '/)
!!!!!!!!!!!!!!!!!!!!!!!!!!

    nout = 6
    if (present(lp)) then
      if (lp < 0) return
      nout = lp
! check lp not equal to the unit number of a file that has been declared
      do i = 1, data%private%nfiles
        if (data%private%unit(i) == lp) then
          nout = 6; exit
        end if
      end do
    end if

    write (nout,'(3a,i3)') ' Error return from OF01_', &
        trim(names(data%entry)), '. Error flag = ', iflag

    if (iflag == -1) then
      write (nout,'(a,i3)') ' Allocation error. stat parameter = ', data%stat
    else if (iflag == -2) then
      write (nout,'(a,i3)') ' Violation of restriction on optional argument '
    else if (iflag == -3) then
      write (nout,'(a,i3)') ' loc out of range'
    else if (iflag == -4) then
      write (nout,'(a,i3)') &
        ' the superfile is not open under OF01'
    else if (iflag == -5) then
      write (nout,'(2a,i3)') ' INQUIRE statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -6) then
      write (nout,'(2a,i3)') ' READ statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -7) then
      write (nout,'(2a,i3)') ' OPEN statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -8) then
      write (nout,'(a,i3)') ' Deallocation error. stat parameter = ', &
        data%stat
    else if (iflag == -9) then
      write (nout,'(a)') &
        ' ifile is out of its range'
    else if (iflag == -10) then
      write (nout,'(2a)') ' one or more superfiles are open through HSL_OF01'
    else if (iflag == -11) then
      write (nout,'(2a)') ' lenw is positive but one or more of the ', &
             'required files does not exist'
    else if (iflag == -12) then
      write (nout,'(2a)') ' filename already exists in path ', &
             trim(data%private%path(data%iostat))
    else if (iflag == -13) then
      write (nout,'(a,i5)') ' file name is longer than ', maxname
    else if (iflag == -14) then
      write (nout,'(2a,i3)') ' CLOSE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -15) then
      write (nout,'(2a,i3)') ' WRITE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -16) then
      write (nout,'(a,i5)') ' path name is longer than ', maxpath
    else if (iflag == -17) then
      write (nout,'(a,i5)') ' unable to open file of given length'
    end if

  end subroutine print_iflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_mcopy(len,buffer,array,map)

! Add buffer into array under the control of map, thus:
!    do l = 1,len
!      i = map(l)
!      array(i) = array(i) + buffer(l)
!    end do

    integer, intent(in) :: len
    real(wp), intent(in) :: buffer(len)
    real(wp), intent(inout) :: array(*)
    integer, intent(in) :: map(len)

! Local variables
    integer :: i,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    k = mod(len,7)
    do l = 1, k
      i = map(l); array(i) = array(i) + buffer(l)
    end do

    do l = k + 1, len, 7
       i = map(l);    array(i) = array(i) + buffer(l)
       i = map(l+1);  array(i) = array(i) + buffer(l+1)
       i = map(l+2);  array(i) = array(i) + buffer(l+2)
       i = map(l+3);  array(i) = array(i) + buffer(l+3)
       i = map(l+4);  array(i) = array(i) + buffer(l+4)
       i = map(l+5);  array(i) = array(i) + buffer(l+5)
       i = map(l+6);  array(i) = array(i) + buffer(l+6)
    end do

  end subroutine of01_mcopy


end module hsl_of01_double

! COPYRIGHT (c) 2007 Council for the Central Laboratory
!               of the Research Councils
! Original date 5 October 2007. Version 1.0.0.
! 13 November 2007. Version 1.1.0 Subroutines of01*copy moved into the module.
!   and removed the separate data module.
! 22 November 2007. Version 2.0.0 Logical argument active replaced by
!    integer(long) argument inactive.
! 7 December 2007. Version 3.0.0 Argument inactive removed from OF01_read
!    and argument retain replaced by discard.
! 29 July 2009. Version 3.1.0 Minor changes made to avoid references to 
!    undefined variables.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hsl_of01_integer
  implicit none
  private

  public of01_data, of01_data_private
  public of01_initialize, of01_open, of01_close, of01_read, of01_write, &
    of01_end


  !!! Parameters !!!
  integer, parameter  :: default_lpage = 2**12 ! Default page length.

  integer, parameter  :: default_maxfiles = 10 ! Maximum number of open files.

  integer, parameter  :: default_npage = 1600 ! Default number of pages
!          in the buffer.

  integer, parameter  :: long = selected_int_kind(18) ! Long integer.

  integer, parameter  :: default_file_size = 2**21 ! Default
!          target length of each file.

  integer, parameter  :: nsup = 2 ! Initial size of the array
!          data%private%filename.

  integer, parameter  :: wp = kind(0) ! Defines data type (single/double).

  integer, parameter  :: ihash = 3 ! Constant used in the hashing function.

  integer, parameter  :: maxpath = 400 ! Max. length of path name.

  integer, parameter  :: maxname = 400 ! Max. length of superfile name.

  integer, parameter  :: maxlen = maxpath+maxname+10 ! Max. length of file name
!          There has to be a limit on the lengths of path and file names because
!          Fortran 95 requires the file name in an OPEN statement to be a scalar
!          character variable.

  integer(wp), parameter :: zero = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type of01_data_private

! Note: Age is measured since last reference.
    integer(wp), allocatable :: buffer(:,:) ! of shape (lpage,npage)

    logical,  allocatable :: differs(:) ! of size npage.
!             differs(I)=.true. if the file version of page I
!             of the buffer is different from the buffer version.

    character(maxname), allocatable :: filename(:)
!             Set to hold a copy of the superfile names given to OF01_open.

    integer(long), allocatable :: first(:) ! of size npage. first(K) is the
!            first page with hash code K or zero if there are none.

    integer :: free !  Start of linked list of free file indices.

    integer(long), allocatable :: left(:),right(:) ! of size maxfiles.
!            A sequence of discarded entries on a superfile is recorded as
!            left(superfile):right(superfile).

    integer(long), allocatable :: highest(:) ! of size maxfiles.
!            highest(I) holds:
!            -1 after call of of01_initialize,
!            -2 after call of of01_close for file I, or
!            >0 highest page number read or written on file I.

    integer, allocatable :: index(:) ! of size npage. index(I) contains the
!            index of the (super)file associated with page I of the buffer.

    integer :: iolength ! iolength of a page

    integer :: maxfiles = 0 ! Maximum number of open files.

    integer, allocatable :: name(:) ! of size maxfiles. For a primary file,
!            name(I) holds the position of the superfile name in filename.

    integer(long), allocatable :: next(:) ! of size npage.
!            next(I) is the next page with same hash code
!            or zero for last in list. Otherwise undefined.

    integer :: nfiles !  Number of file indices in use.

    integer(long) :: nrec ! Number of records in each file.

    integer, allocatable :: nextfile(:) ! of size maxfiles. Used for the
!            linked list of free file indices and linked lists of indices of
!            files in a superfile.  nextfile(I) is the next file index in
!            its list or zero if there are none.

    integer(long), allocatable :: older(:) ! of size npage. older(I) is the
!            next older page to page I, or the youngest page if I is oldest.

    integer(long), allocatable :: page(:) ! of size npage. page(I) is the
!            file address (page number) of page I in the buffer.

    integer(long), allocatable :: page_list(:) ! of size npage.
!            list of pages for this read or write
!            that were in the buffer at the time of call.

    character(maxpath), allocatable :: path(:)
!            The paths given to OF01_initialize

    integer(long), allocatable :: prev(:) ! of size npage. prev(I) is the
!            previous page with same hash code or -(hash code)
!            for first in list. Otherwise undefined.

    integer, allocatable :: unit(:) ! of size maxfiles. unit(I) holds the
!            unit for file I or zero if not in use.

    integer(long), allocatable :: younger(:) ! of size npage.
!            younger(I) is the next younger page to page I, or
!            the oldest page if I is the youngest.

    integer :: youngest ! most recently referenced page in the buffer

  end type of01_data_private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type of01_data

    integer :: entry = 0 ! Index of the entry to the package.

    integer :: iostat    ! Fortran iostat parameter.
    integer :: lpage     ! length of each page of the buffer, that is, the
!            number of scalar variables in each page. The default is 1024.

    integer(long) :: ncall_read  ! number of calls to OF01_read.

    integer(long) :: ncall_write ! number of calls to OF01_write.

    integer(long) :: nio_read ! number of records read by
!            OF01_read and OF01_write.

    integer(long) :: nio_write ! number of records written by
!            OF01_read and OF01_write.

    integer(long) :: npage ! number of pages in the in-core buffer.
!            The default value is 20. It is a long integer, but
!            we do not anticipate it being very large.

    integer(long) :: file_size ! target length of each file.
!             The default value is default_file_size.

    integer(long) :: nwd_read ! number of scalars read by OF01_read.

    integer(long) :: nwd_write ! number of scalars written by OF01_write.

    type (of01_data_private) :: private
! In Fortran 2003, this should be
!   type (of01_data_private), private :: private

    integer :: stat ! Fortran stat parameter.

  end type of01_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  interface of01_initialize
    module procedure of01_initialize_integer
  end interface
  interface of01_open
    module procedure of01_open_integer
  end interface
  interface of01_close
    module procedure of01_close_integer
  end interface
  interface of01_read
    module procedure of01_read_integer
  end interface
  interface of01_write
    module procedure of01_write_integer
  end interface
  interface of01_end
    module procedure of01_end_integer
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_initialize_integer(iflag,data,path,file_size,lpage,npage,lp)

! This subroutine must be called once to initialize
! a structure of derived type of01_data
! and to allocate and initialize its array components.

    integer, intent (out) :: iflag ! Value 0 on successful return.
!            The only possible negative values are:
!            -1 Allocation error. The STAT parameter returned in data%stat.
!            -8 Deallocation error. The STAT parameter returned in data%stat.
!           -16 Character length of path too great.

    type (of01_data) , intent (inout) :: data

    character(*), optional, intent (in) :: path(:) ! Paths for the files.
!            If absent, (\ '' \) is used.

    integer(long), optional, intent (in) :: file_size ! If present and positive,
!            the target length of a file. Default value
!            used if absent or present and not positive.

    integer, optional, intent (in) :: lpage ! If present and positive, the
!            size of each page in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: npage ! If present and positive, the
!            number of pages in the buffer. Default value used
!            if absent or present and not positive.

    integer, optional, intent (in) :: lp ! If present and not negative, the
!            unit number for diagnostic messages. 6 is
!            used if absent or present and negative.

! local variables
    integer :: i !  do loop variable
    integer :: size_path ! size of array path

!!!!!!!!!!!!!!!!!!!!!!!
! Initialise
    iflag = 0
    data%entry       = 1
    data%ncall_read  = 0
    data%ncall_write = 0
    data%nio_read    = 0
    data%nio_write   = 0
    data%nwd_read    = 0
    data%nwd_write   = 0

! If napge/lpage supplied, use instead of the default values
    data%npage = default_npage
    if (present(npage)) then
      if (npage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%npage = npage
    end if

    data%lpage = default_lpage
    if (present(lpage)) then
      if (lpage < 1) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%lpage = lpage
    end if

! Find record length
    if (allocated(data%private%buffer)) then
       deallocate (data%private%buffer,stat=data%stat)
    end if
    allocate (data%private%buffer(data%lpage,1),stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -1; call print_iflag(data,iflag,lp); return
    end if
    data%private%buffer(:,1) = 0
    inquire (iolength=data%private%iolength) data%private%buffer(:,1)
    deallocate (data%private%buffer,stat=data%stat)
    if ( data%stat /= 0) then
       iflag = -8; call print_iflag(data,iflag,lp); return
    end if

    size_path = 1
    if (present(path)) then
       size_path = size(path)
       if (len(path)>maxpath) then
         iflag = -16; call print_iflag(data,iflag,lp); return
       end if
    end if

    data%file_size = default_file_size
    if (present(file_size)) then
      if (file_size < data%lpage) then
        iflag = -2; call print_iflag(data,iflag,lp); return
      end if
      data%file_size = file_size
    end if
    data%private%nrec = data%file_size/data%lpage
    data%file_size = data%private%nrec*data%lpage

! Allocate private arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,   data%private%nextfile, &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%path,      data%private%name,     &
                  data%private%left,      data%private%right,    &
                  data%private%filename, stat=data%stat)
       if (data%stat /= 0) then
         iflag = -8; call print_iflag(data,iflag,lp); return
       end if
    end if
    data%private%maxfiles = default_maxfiles
    allocate (data%private%highest(data%private%maxfiles),  &
              data%private%nextfile(data%private%maxfiles), &
              data%private%index(1:data%npage),             &
              data%private%older(1:data%npage),             &
              data%private%younger(1:data%npage),           &
              data%private%differs(1:data%npage),           &
              data%private%prev(1:data%npage),              &
              data%private%next(1:data%npage),              &
              data%private%first(1:data%npage),             &
              data%private%page(1:data%npage),              &
              data%private%page_list(1:data%npage),         &
              data%private%unit(data%private%maxfiles),     &
              data%private%path(size_path),                 &
              data%private%name(data%private%maxfiles),     &
              data%private%left(data%private%maxfiles),     &
              data%private%right(data%private%maxfiles),    &
              data%private%filename(nsup),                  &
              data%private%buffer(data%lpage,data%npage), stat=data%stat)
    if (data%stat /= 0) then
      data%private%maxfiles = 0
      iflag = -1; call print_iflag(data,iflag,lp); return
    end if

! Initialise arrays
    data%private%free = 0
    data%private%nfiles = 0
    data%private%filename(:) = ''
    do i = 1, data%npage
      data%private%index(i) = -1
      data%private%older(i) = i + 1
      data%private%younger(i) = i - 1
      data%private%page(i) = 0
      data%private%first(i) = 0
      data%private%next(i) = 0
      data%private%prev(i) = 0
      data%private%differs(i) = .false.
    end do
    data%private%youngest = 1
    data%private%younger(1) = data%npage
    data%private%older(data%npage) = 1
    if (present(path)) then
       data%private%path(:) = path(:)
    else
       data%private%path(:) = (/ '' /)
    end if
    data%private%name(:) = 0
    data%private%buffer(:,1) = 0

  end subroutine of01_initialize_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_open_integer(filename,ifile,iflag,data,lenw,lp)

! This subroutine must be called for each superfile that is to be
! accessed through HSL_OF01. It gives the superfile an index and opens
! its files.

    character(len=*), intent (in) :: filename ! Name of the superfile.

    integer, intent (out) :: ifile ! Index of the superfile.

    integer, intent (out) ::  iflag ! Value 0 on successful return.
!             A negative value is associated with an error message
!               on unit lp. Possible negative values are:
!            -1  Allocation error. The STAT parameter returned in data%stat.
!            -5  Error in inquire statement.
!            -7  Error in open statement.
!            -8  Deallocation error. The STAT parameterreturned in data%stat.
!            -11 lenw > 0, but not enough files exist.
!            -12 File filename already exists, but lenw is not present
!                or lenw <=  0.
!            -13 filename is too long.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer(long), optional, intent (in) :: lenw ! length in pages of the
!            part of the file that has been written and is not
!            regarded as having been overwritten by zeros.
!            Pages beyond this are regarded as containing zeros.
!            If not present or lenw <=  0, a new file is opened.

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.
!            If lp < 0, messages are suppressed.

! Local variables
    character(10) :: ci ! Filename extension
    integer :: i ! do loop variable
    integer :: k ! number of secondary files
    integer :: l ! file index
    integer(long) :: lenw_copy ! no. of pages in the virtual array
    integer :: nout ! unit for error messages
    integer :: m ! Previous file index in linked list
    character(maxname), allocatable :: temp(:)

    iflag = 0
    nout = 6
    if (present(lp)) nout = lp
    data%entry = 2

    if (len(filename) > maxname) then
       iflag = -13; go to 100
    end if

! Find the length of the virtual array
    lenw_copy = 0
    if (present(lenw)) lenw_copy = max(0_long,lenw)

! Find number of secondary files
    k = max(0_long,(lenw_copy-1)/data%private%nrec)

! Find suitable indices and units for all the files and open them
    m = 0
    ci=''
    do i = 0,k
      if (lenw_copy > 0) then
        call open_old (data,l,trim(filename)//adjustl(ci),iflag)
      else
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
      end if
      if (iflag /= 0) go to 100
      if (i==0) ifile = l
      data%private%highest(l) = max(0_long,lenw_copy)
      data%private%left(l) = 1
      data%private%right(l) = 0
      if (m > 0) data%private%nextfile(m) = l
      m = l
      data%private%highest(l) = lenw_copy
      write(ci,'(i5)') i + 1
    end do

! Look for a place and store filename there
    k = size(data%private%filename)
    do i = 1, k
      if (data%private%filename(i)== '') then
          data%private%filename(i) = filename
          data%private%name(ifile) = i
          return
      end if
    end do
! Increase size of data%private%filename
    allocate(temp(k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    temp(:) = data%private%filename(:)
    deallocate(data%private%filename,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
    allocate(data%private%filename(2*k),stat=data%stat)
    if (data%stat /= 0) then
      iflag = -1; go to 100
    end if
    data%private%filename(1:k) = temp(:)
    data%private%filename(k+1:) = ''
    deallocate(temp,stat=data%stat)
    if (data%stat /= 0) then
      iflag = -8; go to 100
    end if
! Store filename
    data%private%filename(k+1) = filename
    data%private%name(ifile) = k+1
    return

100 call print_iflag(data,iflag,lp)

  end subroutine of01_open_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_read_integer &
       (ifile,loc,n,read_array,iflag,data,lp,map,discard)

! This subroutine performs the read from a superfile.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be read.

    integer(wp), intent (inout) :: read_array(*)
!            Array into which data from the file is read.
!            If map is present, file data is added into it thus:
!               read_array(map) = read_array(map) + file_data(1:n)

    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not open.
!            -5  Error in Fortran INQUIRE. The IOSTAT parameter is returned in
!                data%iostat.
!            -6  Error in Fortran READ. The IOSTAT parameter is returned in
!                data%iostat.
!            -7  Error in Fortran OPEN statement. The IOSTAT parameter
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

    integer, optional, intent (in) :: map(n) ! map array.

    logical, optional, intent (in) :: discard ! Whether data is to be discarded.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos  ! position within first page of first value
!           required.
    integer(long) :: left,right ! Range of discarded entries
    integer(long) :: ip  ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer
!           pages that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer       :: jfile ! file associated with a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!           required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer       :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: my_discard ! Whether the data is to be discarded.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    if (n <= 0) return
    data%entry = 3

! Check for errors in incoming data
    if (ifile < 1) then
       iflag = -9; go to 200
    else if (ifile > data%private%nfiles) then
       iflag = -4; go to 200
    else if (data%private%highest(ifile) < 0) then
       iflag = -4; go to 200
    else if (loc <= 0) then
       iflag = -3; go to 200
    end if

    my_discard = .false.
    if (present(discard)) then
      if (discard) then
        my_discard = .true.
        left = loc
        right = loc+n-1
        if (loc == data%private%right(ifile)+1) then
          left = data%private%left(ifile)
        else if (loc+n == data%private%left(ifile)) then
            right = data%private%right(ifile)
        end if
        data%private%left(ifile) = left
        data%private%right(ifile) = right
      end if
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_read = data%ncall_read + 1
    data%nwd_read = data%nwd_read + n

    lenl = 0
! Look for required pages that are in the buffer
    look1: do ip = first_page, last_page

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp) == ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
10    call read_from_buffer
      if (my_discard) then
        if ((ip-1)*lpage+1 >= left .and. ip*lpage+1 <= right ) &
          data%private%differs(jp) = .false.
      end if
   end do look1

! look for pages not yet treated and treat each of them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw.le.lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if
! the oldest page is overwritten and becomes the youngest
      jp = data%private%younger(data%private%youngest)
      data%private%youngest = jp

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile>=0) then
! write page from buffer to record page(jp) in the file
! if it has been altered (ie if differs(,jp) is set to 1).
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 200

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip > data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 200
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m > 0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .false.

! perform actual transfer to read_array from buf
      call read_from_buffer

    end do look2
!    call buffer
    return

200 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_from_buffer

      integer :: l1,l2 ! Temporaries
! Read page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip == last_page) l2 = last_pos
      if (.not. present(map)) then
        call FVM_icopy(l2-l1+1,data%private%buffer(l1,jp),read_array(m))
      else
        call of01_mcopy(l2-l1+1,data%private%buffer(l1,jp),read_array,map(m))
      end if
    end subroutine read_from_buffer

  end subroutine of01_read_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_write_integer &
             (ifile,loc,n,write_array,iflag,data,lp,inactive)

! This subroutine performs the writing to a superfile.

    integer, intent (in) :: ifile ! Index of the (super)file.

    integer(long), intent (in) :: loc ! Start position in virtual array.

    integer, intent (in) :: n ! Number of entries to be written.

    integer(wp), intent (in) :: write_array(*) ! array from which data
!            is written to the file.
    integer, intent (out) :: iflag ! Successful return indicated by iflag = 0.
!            Negative value associated with an error message which
!            is output on unit lp. Possible negative values:
!            -3  loc is not positive.
!            -4  Attempt to access a file that is not OF01_open.
!            -5  Error in INQUIRE.  The IOSTAT value
!                is returned in data%iostat.
!            -9  ifile out of range.
!            -15 Error in Fortran WRITE. The IOSTAT parameter
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! unit number for diagnostic messages.
!            Negative for no messages.
!            If not present or equal to the unit number of a file
!            that has already been opened for data, 6 is used.

    integer(long), optional, intent (in) :: inactive ! If present, it is
!            assumed that the data transferred are unlikely to be needed before
!            other data in the buffer. If inactive<loc, it is assumed that this
!            is true also for the data in the interval inactive:loc. If
!            inactive loc+n, it is assumed that this  is true also for the data
!            in the interval loc+n:inactive.

! Local variables
    integer(long) :: first_page ! first page number required.
    integer(long) :: first_pos !  position within first page of first value
!            required.
    integer(long) :: ip ! page number on file of the current page.
    integer(long) :: ipw ! Position in list in page_list(:) of the buffer pages
!            that have been accessed.
    integer(long) :: jh ! temporary variable used to hold a hash code.
    integer(long) :: jp ! page number in the buffer of the required page.
    integer :: jfile ! file index of a page in the buffer
    integer(long) :: l ! temporary variable.
    integer(long) :: last_page ! last page number required.
    integer(long) :: last_pos ! position within last page of last value
!            required.
    integer(long) :: lenl ! length of list in page_list(:) of accessed pages
    integer :: lpage ! page length.
    integer(long) :: m ! temporary variable.
    logical       :: inactive_page ! Whether this page is fully inactive.
    logical       :: inactive_first ! Whether first page is fully inactive.
    logical       :: inactive_middle ! Whether middle pages are fully inactive.
    logical       :: inactive_last ! Whether last page is fully inactive.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    if (n <= 0) return
    data%entry = 4

! Check for errors in incoming data
    if (ifile < 1) then
      iflag = -9; go to 300
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 300
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 300
    else if (loc <= 0) then
      iflag = -3; go to 300
    end if

! lpage is the length of a page (= length of a record in the file)
    lpage = data%lpage

! first_page and last_page are first and last page numbers required.
    first_page = 1 + (loc-1)/lpage
    last_page = 1 + (loc+(n-1)-1)/lpage

    if (present(inactive)) then
       if (first_page==last_page) then
          inactive_first = .false.
       else
          inactive_first = inactive<=1+(first_page-1)*lpage
          inactive_middle = .true.
          inactive_last = inactive>=last_page*lpage
       end if
    else
       inactive_middle = .false.
       inactive_first = .false.
       inactive_last = .false.
    end if

! If there is any overlap with the range of discarded entries,
! reset the range to be null
    if (loc>data%private%right(ifile)) then
    else if (loc+n-1<data%private%left(ifile)) then
    else
       data%private%left(ifile) = 1
       data%private%right(ifile) = 0
    end if

! first_pos and last_pos are positions within first and
! last pages of first and last values required.
    first_pos = loc - (first_page-1)*lpage
    last_pos = loc + n - 1 - (last_page-1)*lpage

    data%ncall_write = data%ncall_write + 1
    data%nwd_write = data%nwd_write + n

    lenl = 0
! Look for required pages in the buffer
    look1: do ip = first_page, last_page
      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if
!     if (inactive_page) write(10,'(i3,i9,3i5)') ifile,loc,n,ip,data%lpage

! See if the page requested is the youngest
      jp = data%private%youngest
      if (data%private%page(jp) == ip) then
        if (data%private%index(jp)==ifile) then
! Special code for when the youngest page is wanted
          lenl = lenl + 1
          data%private%page_list(lenl) = ip
          call write_to_buffer
          go to 10
        end if
      end if
! Find hash code for page ip in file ifile.
      jh = 1 + mod(ip+ifile*ihash,data%npage)
! first(jh) is the first page with hash code jh or 0 if there are none
      jp = data%private%first(jh)
      do l = 1, data%npage+1
        if (jp == 0) cycle look1
! data%private%index(jp) holds the index of the (super)file
! associated with page jp in the buffer
        jfile = data%private%index(jp)
        if (jfile == ifile) then
! data%private%page(jp) is record (page) index in file of page jp in buffer
          if (data%private%page(jp) == ip) exit
        end if
        jp = data%private%next(jp)
      end do
! page ip found as page jp in buffer
      lenl = lenl + 1
      data%private%page_list(lenl) = ip
! remove page from its present position
        l = data%private%younger(jp)
        m = data%private%older(jp)
        data%private%older(l) = m
        data%private%younger(m) = l
! insert page as youngest
        l = data%private%younger(data%private%youngest)
        data%private%older(l) = jp
        data%private%younger(jp) = l
        data%private%older(jp) = data%private%youngest
        data%private%younger(data%private%youngest) = jp
        data%private%youngest = jp
      call write_to_buffer
! If the page is fully inactive, make it the oldest
 10   if (inactive_page) then
         data%private%youngest = data%private%older(jp)
      end if
    end do look1

! look for pages not yet treated and treat them
    ipw = 1
    look2: do ip = first_page, last_page
      if ( ipw <= lenl) then
        if (data%private%page_list(ipw) == ip) then
          ipw = ipw+1
          cycle
        end if
      end if

      if (ip==first_page) then
         inactive_page =  inactive_first
      else if (ip<last_page) then
         inactive_page =  inactive_middle
      else
         inactive_page =  inactive_last
      end if

! the oldest page is overwritten
      jp = data%private%younger(data%private%youngest)
! if transfer is active or the page is partially active, page becomes
! the youngest
     if (.not.inactive_page) then
        data%private%youngest = jp
     end if

! data%private%index(jp) contains index of (super)file associated
! with page jp
      jfile = data%private%index(jp)
      if (jfile >= 0) then
! write page from buffer to record page(jp) in the file
! if it has been altered.
        if (data%private%differs(jp)) call page_write(data, &
          data%private%buffer(:,jp),jfile,data%private%page(jp),iflag)
        if (iflag /= 0) go to 300

! remove page from hash list
! prev(jp) is the last page with same hash code ( or -hash code
! for the first in the list)
        l = data%private%prev(jp)
! next(jp) is next page with same hash code, or 0 if last in list
        m = data%private%next(jp)
        if (m > 0) data%private%prev(m) = l
        if (l > 0) then
          data%private%next(l) = m
        else if (l < 0) then
          l = -l
          data%private%first(l) = m
        end if
      end if

! Read from file record ip into page jp of buffer
! Check whether read extends beyond
! where file written. If it does, fill buffer with zeros and return
      if (ip>data%private%highest(ifile)) then
        data%private%buffer(:,jp) = zero
      else
        call page_read(data,data%private%buffer(:,jp),ifile,ip, &
          iflag)
        if (iflag /= 0) go to 300
      end if
! add to hash list
      jh = 1 + mod(ip+ifile*ihash,data%npage)
      m = data%private%first(jh)
      data%private%next(jp) = m
      if (m>0) data%private%prev(m) = jp
      data%private%first(jh) = jp
      data%private%prev(jp) = -jh
! revise rest of page table
      data%private%index(jp) = ifile
      data%private%page(jp) = ip
      data%private%differs(jp) = .true.

! perform actual transfer to write_array from buf
      call write_to_buffer

    end do look2
!    call buffer
    return

300 call print_iflag(data,iflag,lp)

  contains

   subroutine buffer
     integer i
     i = data%private%youngest
     write(*,'(a)') ' file page'
     do
       if(data%private%index(i)>0) &
       write(*,'(2i5)') data%private%index(i), data%private%page(i)
       i = data%private%older(i)
       if (i == data%private%youngest)exit
     end do
  end  subroutine buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_to_buffer

      integer :: l1,l2 ! Temporaries

! Write page jp of the buffer (page ip of the virtual array) or part thereof.
      l2 = lpage
      l1 = 1
      if (ip == first_page) then
        l1 = first_pos
        m = 1
      else
        m = (ip-first_page)*lpage+2-first_pos
      end if
      if (ip==last_page) l2 = last_pos
      call FVM_icopy(l2-l1+1,write_array(m),data%private%buffer(l1,jp))
      data%private%differs(jp) = .true.

    end subroutine write_to_buffer
  end subroutine of01_write_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_close_integer(ifile,lenw,num_file,iflag,data,lp,lkeep)

! This subroutine should be called when a superfile is finished with.
! Its data in the buffer are unloaded.

    integer, intent (in) :: ifile !  index of the superfile.

    integer(long), intent (out) :: lenw ! length in pages of the part of
!            the file that has been written.
    integer, intent (out) :: num_file ! total number of additional files
!            that have been opened. They are named
!            filenamej, j = 1, ..., num_file.
    integer, intent (out) :: iflag ! iflag = 0 for a successful return.
!            A negative value is associated with an error
!            message on unit lp. Possible negative values:
!            -4  Attempt to access a file that is not OF01_open.
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -9  ifile out of range.
!            -14 Error in CLOSE statement. The IOSTAT value
!                is returned in data%iostat.
!            -15 Error in WRITE.  The IOSTAT value
!                is returned in data%iostat.
!            -17 open unsuccessful for all elements of path.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0 , messages are suppressed.

    logical, optional, intent (in) :: lkeep ! If present and set to .FALSE.,
!            the files are deleted on being closed.
!            Otherwise, the files are kept.

! Local variables
    integer(long) :: i !  do loop variable
    integer :: j       !  do loop variable
    integer :: jfile   !  index of a file
    integer(long) :: l ! page number
    integer :: lpage   ! length of each page (= data%lpage)
    integer(long) :: m ! page number
    integer(long) :: oldest ! oldest page
    integer(long) :: oldest2 ! second oldest page
    character (maxlen) name  ! File name
    character(len=6) :: status ! File status on closing

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 5
! Check for errors
    if (ifile < 1) then
      iflag = -9; go to 400
    else if (ifile > data%private%nfiles) then
      iflag = -4; go to 400
    else if (data%private%highest(ifile) < 0) then
      iflag = -4; go to 400
    end if

    status = 'keep'
    if (present(lkeep)) then
      if (.not. lkeep) status = 'delete'
    end if

! search the pages of the buffer dealing with those associated
! with ifile.
    search:do i = 1, data%npage
! data%private%index(i) is index of the (super)file associated
! with page i of buffer
      if (ifile /= data%private%index(i)) cycle
      data%private%index(i) = -1
      if (data%private%differs(i)) then
! File version of page i is not the same as the buffer version so
! write page i to file
        data%private%differs(i) = .false.
        lpage = data%lpage
        if (status == 'keep') then
          call page_write(data, &
          data%private%buffer(:,i),ifile,data%private%page(i),iflag)
          if (iflag /= 0) then
            call print_iflag(data,iflag,lp)
            return
          end if
        end if
      end if

! remove page from hash list
      l = data%private%prev(i) ! Previous
      m = data%private%next(i) ! Next
      if (m>0) data%private%prev(m) = l
      if (l>0) then
        data%private%next(l) = m
      else if (l < 0) then
        l = -l  ! i is first in the list with code l
        data%private%first(l) = m
      end if

! youngest is the most recently referenced page in the buffer.
! younger(youngest) is the oldest page in the buffer.
      oldest = data%private%younger(data%private%youngest)
! if page i is not already the oldest, make it so.
      if (oldest==i) cycle
! remove page from its present position
      l = data%private%younger(i)
      m = data%private%older(i)
! l is next youngest to page I and m is the next oldest
      data%private%older(l) = m
      data%private%younger(m) = l

! insert page i as the oldest page
! older(oldest) is the youngest page
      data%private%youngest = data%private%older(oldest)

! Page i becomes the oldest page. The page that was the oldest
! becomes the second (ie next) oldest page
      oldest2 = oldest
      oldest = i
      data%private%older(oldest2) = oldest
      data%private%younger(oldest) = oldest2
      data%private%older(oldest) = data%private%youngest
      data%private%younger(data%private%youngest) = oldest
    end do search

    lenw = data%private%highest(ifile)
    if (lenw <= 0 ) status = 'delete'
    num_file = (lenw-1)/data%private%nrec

! remove file from list of declared files.
! close files
    inquire(unit=data%private%unit(ifile),name=name,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -5; go to 400
    end if
    close (unit=data%private%unit(ifile),status=status,iostat=data%iostat)
    if ( data%iostat /= 0) then
      iflag = -14; go to 400
    end if
    data%private%unit(ifile) = 0
    data%private%highest(ifile) = -1
    jfile = ifile
    if (num_file > 0) then
      do j = 1, num_file
        m = data%private%nextfile(jfile)
        inquire(unit=data%private%unit(m),name=name,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -5; go to 400
        end if
        close (unit=data%private%unit(m),status=status,iostat=data%iostat)
        if ( data%iostat /= 0) then
          iflag = -14; go to 400
        end if
        data%private%unit(m) = 0
        data%private%highest(m) = -1
        jfile = m
      end do
    end if
! Put the freed indices at the front of the list list of free indices
    data%private%nextfile(jfile) = data%private%free
    data%private%free = ifile

! Clear the name for reuse
    data%private%filename(data%private%name(ifile))=''

    return

400 call print_iflag(data,iflag,lp)

  end subroutine of01_close_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_end_integer(iflag,data,lp)

! This subroutine deallocates the private array components of data.

    integer, intent (out) :: iflag
!             0  successful return.
!            -8  Error in Fortran DEALLOCATE statement. The STAT
!                value is returned in data%stat.
!            -10 of01_close has not been called for a superfile.

    type (of01_data), intent (inout) :: data

    integer, optional, intent (in) :: lp ! If present, lp must hold the unit
!            for diagnostic messages. Otherwise, lp =  6 used.
!            If lp < 0, messages are suppressed.

! Local variables
    integer :: i

!!!!!!!!!!!!!!!!!!!
    iflag = 0
    data%entry = 6

    do i = 1, data%private%nfiles
      if (data%private%highest(i) >= 0) then
        iflag = -10; call print_iflag(data,iflag,lp)
        return
      end if
    end do

! Deallocate arrays
    if ( allocated(data%private%highest)) then
      deallocate (data%private%highest,  data%private%nextfile,  &
                  data%private%index,     data%private%older,    &
                  data%private%younger,   data%private%differs,  &
                  data%private%prev,      data%private%next,     &
                  data%private%first,     data%private%page,     &
                  data%private%page_list, data%private%unit,     &
                  data%private%buffer,    data%private%path,     &
                  data%private%left,      data%private%right,    &
                  data%private%name,      data%private%filename, &
                  stat=data%stat)
      data%private%maxfiles = 0
      if (data%stat /= 0) then
        iflag = -8; call print_iflag(data,iflag,lp)
      end if
    end if

  end subroutine of01_end_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_write(data,buffer,ifile,locp,iflag)

! This subroutine transfers a page of the buffer to a file. If the
! page is within the existing virtual array then a simple write
! is performed. To write a page beyond the present end may require
! writing pages full of zeros in front of its position.

    type (of01_data), intent (inout) :: data

    integer(wp), intent (inout) :: buffer(data%lpage) ! array from which
!             the transfer is made.
    integer, intent (in) :: ifile ! the index of the primary file.

    integer(long), intent (in) :: locp ! page number in the file.

    integer, intent (out) :: iflag ! error flag. Possible negative values:
!            -5  Error in INQUIRE. IOSTAT value is returned in data%iostat.
!            -7  Error in OPEN. IOSTAT valueis returned in data%iostat.
!            -15 Error in WRITE. IOSTAT value is returned in data%iostat.

! Local variables
    integer(long) :: hpage ! highest page read or written for ifile
    integer(long) :: l ! do loop variable

    hpage = data%private%highest(ifile)
    if (locp <= hpage+1) then
! Write straight to file
      call actual_page_write(locp)
      if (iflag /= 0) return
    else
! Fill the gap with zero page
! first write out buffer to page hpage+1
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
      buffer(:) = zero
! write pages of zeros
      do l = hpage+2, locp-1
        call actual_page_write(l)
        if (iflag /= 0) return
      end do
! read buffer back in from page hpage+1
      call page_read(data,buffer,ifile,hpage+1,iflag)
      if (iflag /= 0) return
! Now write buffer to the required location
      call actual_page_write(locp)
      if (iflag /= 0) return
! Finally write zeros to hpage+1
      buffer(:) = zero
      call actual_page_write(hpage+1)
      if (iflag /= 0) return
    end if

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine actual_page_write(locp)

! Perform an actual i/o write of a page to a superfile. If necessary, open
! intermediate secondary files.

    integer(long), intent (in) :: locp ! page index in the superfile.

! Local variables
    character(10) :: ci ! Filename extension.
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or >0 for k-th secondary file.
    integer l ! index of the file that is written.
    integer m ! Previous file in linked list.
    character(len=maxlen) :: filename ! name of the primary file.

    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec
    l = ifile
    do i = 1, k
      m = l
      l = data%private%nextfile(l)
      if (l == 0) then
! open new file (filename1, filename2, etc)
        filename = data%private%filename(data%private%name(ifile))
        write(ci,'(i5)')i
        call open_new (data,l,trim(filename)//adjustl(ci),iflag)
        if (iflag /= 0)  return

! Add to linked list for ifile
        if (m > 0) data%private%nextfile(m) = l
        m = l
      end if
    end do

    data%nio_write = data%nio_write + 1
    write (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
!    write(*,*) 'write record',new_locp,data%iostat
    if (data%iostat /= 0) then
      iflag = -15; return
    end if
    data%private%highest(ifile) = max(data%private%highest(ifile),locp)

  end subroutine actual_page_write

  end subroutine page_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine page_read(data,buffer,ifile,locp,iflag)

! This subroutine performs an actual i/o read of a page from a superfile

    type (of01_data), intent (inout) :: data

    integer(wp), intent (out) :: buffer(data%lpage) ! array to receive page.

    integer, intent (in) :: ifile ! index of the superfile.

    integer(long), intent (in) :: locp ! page index in the superfile.

    integer, intent (out) :: iflag ! Error flag. Possible negative value:
!            -6 error in read statement

! Local variables
    integer(long) :: new_locp ! page index in the file that is read.
    integer i ! Do index.
    integer k ! 0 for primary file or > 0 for k-th secondary file.
    integer l ! index of the file that is read.

!!!!!!!!!!!!!!!!!!!!!!!!!!
    iflag = 0
    k = (locp-1)/data%private%nrec
    new_locp = locp - k*data%private%nrec

! Find index of the kth linked file
    l = ifile
    do i = 1, k
      l = data%private%nextfile(l)
    end do

    data%nio_read = data%nio_read + 1
    read (unit=data%private%unit(l),rec=new_locp,iostat=data%iostat) buffer
    if (data%iostat /= 0) iflag = -6

  end subroutine page_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_new (data,ifile,filename,iflag)

! Open a new file

    type (of01_data) :: data

    integer, intent(out) :: ifile ! File index.

    character (*) :: filename!    ! File name.

    integer, intent (inout) :: iflag ! Error flag. Possible nonzero value:
!            -7 Error in OPEN.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: j ! Do index.
    integer :: unit ! Unit.
    logical :: open ! open variable for inquire.

! Check file does not already exist
    do i = 1, size(data%private%path)
       inquire (file=trim(data%private%path(i))//filename,exist=ex,&
          iostat=data%iostat)
       if (data%iostat /= 0) then
          iflag = -5; return
       else if (ex) then
          iflag = -12; data%iostat=i; return
       end if
    end do

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

! Open the file on the first available path
    outer: do i = 1, size(data%private%path)
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       if (data%iostat /= 0) then
         close(unit,iostat=data%iostat,status='delete')
         cycle
       end if

       if (size(data%private%path) == 1) return
! Check that there is room for the file
       do j = 1, data%private%nrec
         write(unit,rec=j,iostat=data%iostat)data%private%buffer(:,1)
         if (data%iostat /= 0) then
           close(unit,iostat=data%iostat,status='delete')
           cycle outer
         end if
       end do
       close(unit,iostat=data%iostat)
       if (data%iostat /= 0) then
         inquire(unit,opened=open)
         if (.not.open) &
         open (unit,file=trim(data%private%path(i))//filename, &
             access='direct',iostat=data%iostat,recl=data%private%iolength)
         close(unit,iostat=data%iostat,status='delete')
         cycle outer
       end if
       open (unit,file=trim(data%private%path(i))//filename, &
          access='direct',iostat=data%iostat,recl=data%private%iolength)
       return
    end do outer
    iflag = -17

  end subroutine open_new

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine open_old (data,ifile,filename,iflag)

! Open an old file.

    type (of01_data) :: data

    integer, intent(out) :: ifile ! File index.
    character (*) :: filename     ! File name.

    integer, intent (inout) :: iflag ! Possible nonzero values:
!            -5  Error in inquire statement
!            -7  Error in OPEN.  The IOSTAT value is returned in data%iostat.
!            -11 File does not exist

! Local variables
    logical :: ex ! exist variable for inquire.
    integer :: i ! Do index.
    integer :: unit ! Unit.

! Find the file.
    do i = 1, size(data%private%path)
      inquire (file=trim(data%private%path(i))//filename,exist=ex,&
             iostat=data%iostat)
      if (data%stat /= 0) then
        iflag = -5; return
      end if
      if (ex) exit
    end do
    if (.not. ex) then
! File not found
      iflag = -11; return
    end if

    call find_unit(data,unit,ifile,iflag)
    if (iflag /= 0) return

    open (unit,file=trim(data%private%path(i))//filename, access='direct', &
          iostat=data%iostat, recl=data%private%iolength)
    if (data%iostat /= 0) then
      iflag = -7; return
    end if

  end subroutine open_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_unit(data,unit,ifile,iflag)

! Find a free unit and a free file index.

    type (of01_data) :: data
    integer, intent (out) :: unit ! Free unit.

    integer, intent (out) :: ifile ! Free file.

    integer, intent (inout) :: iflag ! Possible nonzero value:
!           -5 Error in INQUIRE.  The IOSTAT value is returned in data%iostat.

! Local variables
    logical :: open
    integer :: oldsize

! Find free unit
    do unit = 12, huge(0)
      inquire (unit=unit,iostat=data%iostat,opened=open)
      if (data%iostat /= 0) then
! Error in inquire
        ifile = -1
        iflag = -5; return
      end if
      if (.not. open) exit
    end do

! Find free file index
    ifile = data%private%free
    if (ifile > 0) then
      data%private%free = data%private%nextfile(ifile)
    else
      ifile =  data%private%nfiles+1
      if (ifile>data%private%maxfiles) then
! Reallocate files
        data%private%maxfiles = 0
        call reallocate_long(data%private%highest)
        if (iflag < 0) return
        call reallocate_long(data%private%left)
        if (iflag < 0) return
        call reallocate_long(data%private%right)
        if (iflag < 0) return
        call reallocate(data%private%nextfile)
        if (iflag < 0) return
        oldsize = size(data%private%name)
        call reallocate(data%private%name)
        data%private%name(oldsize+1:) = 0
        if (iflag < 0) return
        call reallocate(data%private%unit)
        if (iflag < 0) return
        data%private%maxfiles = size(data%private%unit)
      end if
      data%private%nfiles =  ifile
    end if
    data%private%nextfile(ifile) = 0
    data%private%highest(ifile) = -1
    data%private%left(ifile) = 1
    data%private%right(ifile) = 0
    data%private%unit(ifile) = unit

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    subroutine reallocate(array)

      integer, allocatable :: array(:)

      integer, allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)
      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine reallocate_long(array)

      integer(long), allocatable :: array(:)

      integer(long), allocatable :: temp(:)

      allocate(temp(size(array)),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      temp(:) = array(:)

      deallocate(array,stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -8; return
      end if

      allocate(array(size(temp)*2),stat=data%stat)
      if ( data%stat /= 0) then
        iflag = -1; return
      end if

      array(1:size(temp)) = temp(:)

      deallocate(temp,stat=data%stat)
      if ( data%stat /= 0) then
       iflag = -8; return
      end if

    end subroutine reallocate_long

  end subroutine find_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine print_iflag(data,iflag,lp)

    type (of01_data) :: data

    integer, intent (in) :: iflag ! Error flag.

    integer, optional :: lp ! unit number for diagnostic messages.
!            Negative for no messages. If not present or equal to unit
!            number of a file that has already been opened for data, 6 is used.

! Local variables
    integer :: i, nout
    character(10), parameter :: names(6) = (/'initialize','open      ',&
                    'read      ','write     ','close     ','end       '/)
!!!!!!!!!!!!!!!!!!!!!!!!!!

    nout = 6
    if (present(lp)) then
      if (lp < 0) return
      nout = lp
! check lp not equal to the unit number of a file that has been declared
      do i = 1, data%private%nfiles
        if (data%private%unit(i) == lp) then
          nout = 6; exit
        end if
      end do
    end if

    write (nout,'(3a,i3)') ' Error return from OF01_', &
        trim(names(data%entry)), '. Error flag = ', iflag

    if (iflag == -1) then
      write (nout,'(a,i3)') ' Allocation error. stat parameter = ', data%stat
    else if (iflag == -2) then
      write (nout,'(a,i3)') ' Violation of restriction on optional argument '
    else if (iflag == -3) then
      write (nout,'(a,i3)') ' loc out of range'
    else if (iflag == -4) then
      write (nout,'(a,i3)') &
        ' the superfile is not open under OF01'
    else if (iflag == -5) then
      write (nout,'(2a,i3)') ' INQUIRE statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -6) then
      write (nout,'(2a,i3)') ' READ statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -7) then
      write (nout,'(2a,i3)') ' OPEN statement error;', &
                             ' iostat parameter = ', data%iostat
    else if (iflag == -8) then
      write (nout,'(a,i3)') ' Deallocation error. stat parameter = ', &
        data%stat
    else if (iflag == -9) then
      write (nout,'(a)') &
        ' ifile is out of its range'
    else if (iflag == -10) then
      write (nout,'(2a)') ' one or more superfiles are open through HSL_OF01'
    else if (iflag == -11) then
      write (nout,'(2a)') ' lenw is positive but one or more of the ', &
             'required files does not exist'
    else if (iflag == -12) then
      write (nout,'(2a)') ' filename already exists in path ', &
             trim(data%private%path(data%iostat))
    else if (iflag == -13) then
      write (nout,'(a,i5)') ' file name is longer than ', maxname
    else if (iflag == -14) then
      write (nout,'(2a,i3)') ' CLOSE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -15) then
      write (nout,'(2a,i3)') ' WRITE statement error;', &
                            ' iostat parameter = ', data%iostat
    else if (iflag == -16) then
      write (nout,'(a,i5)') ' path name is longer than ', maxpath
    else if (iflag == -17) then
      write (nout,'(a,i5)') ' unable to open file of given length'
    end if

  end subroutine print_iflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine FVM_icopy(len,buffer,array)
! Copy buffer into array, thus:
!    do l = 1,len
!      array(i) =  buffer(i)
!    end do

    integer, intent(in) :: len
    integer(wp), intent(in) ::  buffer(len)
    integer(wp), intent(inout) :: array(*)

    integer :: i,k
    k = mod(len,7)
    do i = 1, k
      array(i) =  buffer(i)
    end do

    do i = k + 1, len, 7
       array(i)   = buffer(i)
       array(i+1) = buffer(i+1)
       array(i+2) = buffer(i+2)
       array(i+3) = buffer(i+3)
       array(i+4) = buffer(i+4)
       array(i+5) = buffer(i+5)
       array(i+6) = buffer(i+6)
    end do

  end subroutine FVM_icopy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine of01_mcopy(len,buffer,array,map)

! Add buffer into array under the control of map, thus:
!    do l = 1,len
!      i = map(l)
!      array(i) = array(i) + buffer(l)
!    end do

    integer, intent(in) :: len
    integer(wp), intent(in) :: buffer(len)
    integer(wp), intent(inout) :: array(*)
    integer, intent(in) :: map(len)

! Local variables
    integer :: i,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    k = mod(len,7)
    do l = 1, k
      i = map(l); array(i) = array(i) + buffer(l)
    end do

    do l = k + 1, len, 7
       i = map(l);    array(i) = array(i) + buffer(l)
       i = map(l+1);  array(i) = array(i) + buffer(l+1)
       i = map(l+2);  array(i) = array(i) + buffer(l+2)
       i = map(l+3);  array(i) = array(i) + buffer(l+3)
       i = map(l+4);  array(i) = array(i) + buffer(l+4)
       i = map(l+5);  array(i) = array(i) + buffer(l+5)
       i = map(l+6);  array(i) = array(i) + buffer(l+6)
    end do

  end subroutine of01_mcopy


end module hsl_of01_integer

! COPYRIGHT (c) 2007 Science and Technology Facilities Council
! Original date 21st Sept 2007

!-*-*-*-*-*-  HSL_KB22_ long_integer MODULE -*-*-*-*-*-

!  Nick Gould and Philippe Toint, for GALAHAD productions
!  Copyright reserved
!  January 26th 1995

! 21st September 2007 Version 1.0.0. Version numbering added.

      MODULE HSL_KB22_long_integer

         IMPLICIT NONE

         PRIVATE
         PUBLIC :: KB22_build_heap, KB22_get_smallest

!-----------------------------
!   L o n g   i n t e g e r
!-----------------------------

         INTEGER, PARAMETER :: long = selected_int_kind( 18 )

         INTERFACE KB22_build_heap
             MODULE PROCEDURE KB22_build_heap_long_integer
         END INTERFACE
         INTERFACE KB22_get_smallest
             MODULE PROCEDURE KB22_get_smallest_long_integer
         END INTERFACE

      CONTAINS

!-*-*-*-  H S L _ K B 2 2 _ b u i l d _ h e a p   S U B R O U T I N E   -*-*-*

         SUBROUTINE KB22_build_heap_long_integer( n, A, inform, INDA )

!  Given an array A, elements A(1), ...., A(N), subroutine KB22_build_heap
!  re-arranges the elements to form a heap in which each parent has a smaller
!  value than either of its children.

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures SETHEAP and INHEAP

!  ------------------------- dummy arguments --------------------------
!
!  n      integer, which gives the number of values to be sorted.
!         n must be positive
!
!  A      long integer array of length n. On input, A must contain the
!         values which are to be sorted. On output, these values
!         will have been permuted so as to form a heap
!
!  inform integer, which informs the user of the success of KB22_build_heap.
!         If inform = 0 on exit, the heap has been formed.
!         If inform = 1 on exit, n was input with a value less than
!                       or equal to 0 and the heap has not been formed.
!
!  INDA   is an OPTIONAL integer array of length n. On input, INDA may
!         be used to hold indexing information (such as the original
!         order of the values) about A. On output, INDA will have been
!         permuted so that INDA(k) still refers to A(k).
!
!  ------------------ end of dummy arguments --------------------------

         INTEGER, INTENT( IN ) :: n
         INTEGER, INTENT( OUT ) :: inform
         INTEGER ( KIND = long ), INTENT( INOUT ), DIMENSION( n ) :: A
         INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( n ) :: INDA

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: i, j, k, indin
         INTEGER ( KIND = long ) :: rin
         LOGICAL :: index

!  Add the elements to the heap one at a time

         index = PRESENT( INDA )
         IF ( n <= 0 ) THEN
            inform = 1
            RETURN
         ENDIF

         IF ( index ) THEN
           DO k = 2, n
              rin = A( k )
              indin = INDA( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

              i = k
              DO
                 IF ( i <= 1 ) EXIT
                 j = i / 2
                 IF ( A( j ) <= rin ) EXIT
                 A( i ) = A( j )
                 INDA( i ) = INDA( j )
                 i = j
              END DO
              A( i ) = rin
              INDA( i ) = indin
           END DO
         ELSE
           DO k = 2, n
              rin = A( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

              i = k
              DO
                 IF ( i <= 1 ) EXIT
                 j = i / 2
                 IF ( A( j ) <= rin ) EXIT
                 A( i ) = A( j )
                 i = j
              END DO
              A( i ) = rin
           END DO
         END IF
         inform = 0

         RETURN

!  End of subroutine KB22_build_heap

         END SUBROUTINE KB22_build_heap_long_integer

!-*-*-  H S L _ K B 2 2 _ g e t _ s m a l l e s t  S U B R O U T I N E   -*-*-

         SUBROUTINE KB22_get_smallest_long_integer( m, A, inform, INDA )

!  Given an array A, elements A(1), ...., A(m) forming a heap,
!  KB22_get_smallest assigns to rout the value of A(1), the smallest
!  member of the heap, and arranges the remaining members as elements
!  1 to m - 1 of A. rout is then placed in A(m)

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures OUTHEAP and SWOPHEAP

!  ------------------------- dummy arguments --------------------------
!
!  m      integer, which gives the number of values to be sorted.
!         m must be positive
!
!  A      long integer array of length m. On input, A must contain the values
!         which are to be sorted stored in a heap. On output, the smallest
!         value will have been moved into A(m) and the remaining values A(k),
!         k = 1,..., m-1 will have been restored to a heap
!
!  inform integer, which informs the user of the success of KB22_get_smallest.
!         If inform = 0 on exit, the smallest value has been found.
!         If inform = 1 on exit, m was input with a value less than
!                       or equal to 0 and the heap has not been formed
!
!  INDA   optional integer array of length m. On input, INDA may be used
!         to hold indexing information (see KB22_build_heap) about A.
!         On output, INDA will have been permuted so that INDA(k) still
!         refers to A(k). This argument is only permitted if it was present
!         when calling KB22_build_heap
!
!  ------------------ end of dummy arguments --------------------------

         INTEGER, INTENT( IN ) :: m
         INTEGER, INTENT( OUT ) :: inform
         INTEGER ( KIND = long ), INTENT( INOUT ), DIMENSION( m ) :: A
         INTEGER, INTENT( INOUT ), OPTIONAL, DIMENSION( m ) :: INDA

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

         INTEGER :: i, j, indin, indout
         INTEGER ( KIND = long ) :: rin, rout
         LOGICAL :: index

         index = PRESENT( INDA )

!  Add the element rin to the heap, extract and assign to rout
!  the value of the smallest member of the resulting set, and
!  leave the remaining elements in a heap of the original size.
!  In this process, elements 1 to n+1 of the array A may be disturbed

         IF ( m <= 0 ) THEN
            inform = 1
            RETURN
         ENDIF

         IF ( m > 1 ) THEN

           IF ( index ) THEN
              i = 1
              rout = A( 1 )
              indout = INDA( 1 )
              rin = A( m )
              indin = INDA( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

              DO
                 j = i + i
                 IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

                 IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

                 IF ( A( j ) >= rin ) EXIT
                 A( i ) = A( j )
                 INDA( i ) = INDA( j )
                 i = j
              END DO

!  The heap has been restored

              A( i ) = rin
              INDA( i ) = indin

!  Store the smallest value in the now vacated m-th position of the list

              A( m ) = rout
              INDA( m ) = indout
           ELSE
              i = 1
              rout = A( 1 )
              rin = A( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

              DO
                 j = i + i
                 IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

                 IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

                 IF ( A( j ) >= rin ) EXIT
                 A( i ) = A( j )
                 i = j
              END DO

!  The heap has been restored

              A( i ) = rin

!  Store the smallest value in the now vacated m-th position of the list

              A( m ) = rout
           END IF

         END IF
         inform = 0

         RETURN

!  End of subroutine KB22_get_smallest

         END SUBROUTINE KB22_get_smallest_long_integer

!  End of module HSL_KB22_long_integer

      END MODULE HSL_KB22_long_integer
! (c) STFC 2010--2011
! Originating author: Jonathan Hogg
!
! Given a pivot order, this package performs common tasks
! required in the analyse phase of a symmetric sparse direct solver.
! Either the entire analyse may be performed or individual tasks.
! The matrix may be hled in assembled form or in elemental form.
!
! Version 1.1.0
! See ChangeLog for version history

! To convert to long:
! s/_long/_long
! Set pkg_type to long
module hsl_mc78_integer
   implicit none

   private
   public :: mc78_control
   public :: mc78_analyse, mc78_supervars, mc78_compress_by_svar, mc78_etree, &
      mc78_elt_equiv_etree, mc78_postorder, mc78_col_counts, mc78_supernodes, &
      mc78_stats, mc78_row_lists, mc78_optimize_locality

   integer, parameter :: dp = kind(0d0) ! not package type
   integer, parameter :: long = selected_int_kind(18)

   integer, parameter :: minsz_ms = 16 ! minimum size to use merge sort

   integer, parameter :: pkg_type = kind(0) ! package type - integer or long

   type mc78_control
      integer :: heuristic = 1 ! 1=ma77 2=cholmod
      integer :: nrelax(3) = (/ 4, 16, 48 /) ! CHOLMOD-like
      real(dp) :: zrelax(3) = (/ 0.8, 0.1, 0.05 /) ! CHOLMOD-like
      integer :: nemin = 16  ! Node amalgamation parameter

      integer :: unit_error = 6
      integer :: unit_warning = 6
      logical :: ssa_abort = .false. ! If .true., then return with an error if
         ! an assembled matrix is detected as symbolically singular (we do
         ! not garuntee to detect _all_ symbolically singular matrices).
         ! If .false., then a warning is raised instead.

      logical :: svar = .false. ! If .true. then supervariables are used in
         ! the assembled case, otherwise they are not. Supervaraibles are
         ! always used in the elemental case.
      logical :: sort = .false. ! If .true. then entries within each supernode's
         ! row lists are sorted. Otherwise they might not be.
      logical :: lopt = .false. ! If .true. then variable ordering is optimized
         ! for cache locality. Otherwise it is not.
   end type mc78_control

   integer, parameter :: MC78_ERROR_ALLOC = -1 ! allocate error
   integer, parameter :: MC78_ERROR_SSA   = -2 ! symbolically singular assembled
   integer, parameter :: MC78_ERROR_ROW_SMALL = -3 ! supplied row array to short
   integer, parameter :: MC78_ERROR_UNKNOWN = -99 ! internal/unknown error

   ! Warning flags are treated as bit masks, add together if multiple occour
   integer, parameter :: MC78_WARNING_SSA = 1 ! symbolically singular assembled
   integer, parameter :: MC78_WARNING_BLK_SVAR = 2 ! svar and blk pivs requested

   interface mc78_analyse
      module procedure mc78_analyse_assembled_integer
      module procedure mc78_analyse_elemental_integer
   end interface mc78_analyse

   interface mc78_supervars
      module procedure mc78_supervars_integer
   end interface mc78_supervars

   interface mc78_compress_by_svar
      module procedure mc78_compress_by_svar_integer
   end interface mc78_compress_by_svar

   interface mc78_etree
      module procedure mc78_etree_integer
   end interface mc78_etree

   interface mc78_elt_equiv_etree
      module procedure mc78_elt_equiv_etree_integer
   end interface

   interface mc78_postorder
      ! Note: cannot distinguish postorder_std between integer and long versions
      module procedure mc78_postorder_std
      module procedure mc78_postorder_detect
   end interface mc78_postorder

   interface mc78_col_counts
      module procedure mc78_col_counts_integer
   end interface mc78_col_counts

   ! Note: cannot distinguish mc78_supernodes between integer and long versions
   ! Note: cannot distinguish mc78_stats between integer and long versions

   interface mc78_row_lists
      module procedure mc78_row_lists_nosvar_integer
      module procedure mc78_row_lists_svar_integer
   end interface mc78_row_lists

   ! Note: cannot distinguish mc78_optimize_locality between integer and long
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Main analysis routines   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! For assembled matrix input, this subroutine performs a full analysis.
! This is essentially a wrapper around the rest of the package.
!
! Performance might be improved by:
! * Improving the sort algorithm used in find_row_idx
!
subroutine mc78_analyse_assembled_integer(n, ptr, row, perm, nnodes, sptr, &
      sparent, rptr, rlist, control, info, stat, nfact, nflops, piv_size)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer(pkg_type), dimension(:), allocatable :: bptr ! copy of matrix with
      ! added entries for block pivots - column pointers
   integer, dimension(:), allocatable :: brow ! copy of matrix with added
      ! entries for block pivots - row indices
   integer :: flag ! return status flag for call to compress_by_svar
   integer(pkg_type) :: sz
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer, dimension(:), allocatable :: sinvp
   integer :: j
   integer :: k
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: sperm
   integer(pkg_type), dimension(:), allocatable :: ptr2
   integer :: realn ! number of variables with an actual entry present
   integer, dimension(:), allocatable :: row2
   integer :: st ! stat argument in allocate calls
   logical :: svar_r

   integer :: svar_type ! 0=none, 1=col, 2=compressed form

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   svar_r = control%svar

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) then
      call convert_to_blk_piv(n, invp, piv_size)
      allocate(bptr(n+1), brow(ptr(n+1)-1+2*n), stat=st)
      if(st.ne.0) goto 490
      call mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, piv_size, &
         st)
      if(st.ne.0) goto 490
      if(svar_r) then
         ! Supervariables don't interact well with block pivots, so don't do it
         svar_r = .false.
         info = info + MC78_WARNING_BLK_SVAR
      endif
   endif

   ! Determine supervariables (if required)
   if(svar_r) then
      allocate(svara(n), stat=st)
      if(st.ne.0) goto 490
      realn = n
      call mc78_supervars(realn, ptr, row, perm, invp, nsvar, svara, st)
      if(st.ne.0) goto 490
      if(n.ne.realn) then
         if(control%ssa_abort) then
            if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
               "HSL_MC78: Error, matrix is symbolically singular and ", &
               "control%ssa_abort=.true.."
            info = MC78_ERROR_SSA
            return
         else
            if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
               "HSL_MC78: Warning, matrix is symbolically singular."
            info = info + MC78_WARNING_SSA
         endif
      endif
      if(3*nsvar.lt.2*n) then
         svar_type = 2 ! Use compressed form
      else
         svar_type = 0 ! Do not use supervariables
      endif

      select case(svar_type)
      case(0) ! do not use supervariables
         ! release resources
         deallocate(svara, stat=st)
      case(2) ! Compressed form
         ! It is worth using the compressed form
         ! Determine upper bound on size of data for compressed array
         sz = 0
         k = 1
         do i = 1, nsvar
            j = invp(k)
            sz = sz + ptr(j+1) - ptr(j)
            k = k + svara(i)
         end do
         allocate(ptr2(nsvar+1), row2(sz), sperm(nsvar), sinvp(nsvar), stat=st)
         if(st.ne.0) goto 490
         call mc78_compress_by_svar(n, ptr, row, invp, nsvar, svara, &
            ptr2, sz, row2, flag, st)
         select case(flag)
         case(0) ! Everything OK
            ! Do nothing
         case(-1) ! Allocate failure
            goto 490
         case default ! Should never happen
            info = MC78_ERROR_UNKNOWN
            return
         end select
         ! Compressed matrix is in pivot order already
         do i = 1, nsvar
            sperm(i) = i
            sinvp(i) = i
         end do
      end select
   else
      svar_type = 0
      realn = n ! Assume full rank
   endif

   select case(svar_type)
   case(0)
      if(present(piv_size)) then
         call mc78_inner_analyse(n, realn, bptr, brow, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st, &
            block_pivots=piv_size)
      else
         call mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, &
            sptr, sparent, scc, rptr, rlist, control, info, st)
      endif
   case(2)
      call mc78_inner_analyse(nsvar, i, ptr2, row2, sperm, sinvp, nnodes, &
         sptr, sparent, scc, rptr, rlist, control, info, st, wt=svara, &
         block_pivots=piv_size)
      if(st.ne.0) goto 490
      if(info.lt.0) return
      if(i.ne.nsvar) then
         ! Note: This code should NEVER execute
         if(control%unit_error.gt.0) &
            write(control%unit_error, "(a,2(a,i8))") "MC78_ANALYSE Internal ", &
               "Error: supervariable matrix is rank deficient: i = ", i, &
               "nsvar = ", nsvar
         info = MC78_ERROR_UNKNOWN
         return
      endif
      call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, sptr, st)
   end select
   if(st.ne.0) goto 490
   if(info.lt.0) return

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn = ", n, realn
   !print *, "ptr = ", ptr
   !do i = 1, n
   !   print *, "row(", i, ") = ", row(ptr(i):ptr(i+1)-1)
   !end do
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm
   !print *, "piv_size = ", piv_size

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_assembled_integer

!
! Inner core for assembled analyse routine, used to make calls in compressed
! (supervariable) case and standard case uniform
!
subroutine mc78_inner_analyse(n, realn, ptr, row, perm, invp, nnodes, sptr, &
      sparent, scc, rptr, rlist, control, info, st, wt, block_pivots)
   integer, intent(in) :: n ! Dimension of system
   integer, intent(out) :: realn ! Symbolic dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer, dimension(:), allocatable, intent(out) :: scc ! supernodal col cnt
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st
   integer, dimension(n), optional, intent(in) :: wt ! Weights of columns
      ! (i.e. size of each supervariable they represent)
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: ntot ! total number of variables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: tperm ! temporary permutation vector

   ! Build elimination tree
   allocate(parent(n), stat=st)
   if(st.ne.0) return
   call mc78_etree(n, ptr, row, perm, invp, parent, st)
   if(st.ne.0) return

   ! Postorder tree (modifies perm!)
   call mc78_postorder(n, realn, ptr, perm, invp, parent, st, block_pivots)
   if(st.ne.0) return

   if(n.ne.realn) then
      if(control%ssa_abort) then
         if(control%unit_error.gt.0) write(control%unit_error, "(2a)") &
            "HSL_MC78: Error, matrix is symbolically singular and ", &
            "control%ssa_abort=.true.."
         info = MC78_ERROR_SSA
         return
      else
         if(control%unit_warning.gt.0) write(control%unit_warning, "(a)") &
            "HSL_MC78: Warning, matrix is symbolically singular."
         info = info + MC78_WARNING_SSA
      endif
   endif

   ! Determine column counts
   allocate(cc(n+1), stat=st)
   if(st.ne.0) return
   call mc78_col_counts(n, ptr, row, perm, invp, parent, cc, st, wt=wt)
   if(st.ne.0) return

   ! Identify supernodes
   allocate(tperm(n), sptr(n+1), sparent(n), scc(n), stat=st)
   if(st.ne.0) return
   call mc78_supernodes(n, realn, parent, cc, tperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt=wt, block_pivots=block_pivots)
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(n, tperm, perm, invp, cc, block_pivots=block_pivots)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum(scc(1:nnodes))), stat=st)
   if(st.ne.0) return
   if(present(wt)) then
      ntot = sum(wt)
      call mc78_row_lists(n, wt, ntot, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   else
      call mc78_row_lists(n, ptr, row, perm, invp, nnodes, sptr, &
         sparent, scc, rptr, rlist, control, info, st)
   endif
   if(st.ne.0) return
end subroutine mc78_inner_analyse

!
! This subroutine performs full analyse when A is in elemental form.
! This is essentially a wrapper around the rest of the package.
!
subroutine mc78_analyse_elemental_integer(n, nelt, starts, vars, perm, &
      eparent, nnodes, sptr, sparent, rptr, rlist, control, info, stat, &
      nfact, nflops, piv_size)
   integer, intent(in) :: n ! Maximum integer used to index an element
   integer, intent(in) :: nelt ! Number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! Element pointers
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! Variables
      !assoicated with each element. Element i has vars(starts(i):starts(i+1)-1)
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(nelt), intent(out) :: eparent ! On exit, eparent(i) holds
      ! node of assembly that element i is a child of.
   integer, intent(out) :: nnodes ! number of supernodes found
   integer, dimension(:), allocatable, intent(out) :: sptr ! supernode pointers
   integer, dimension(:), allocatable, intent(out) :: sparent ! assembly tree
   integer(long), dimension(:), allocatable, intent(out) :: rptr
      ! pointers to rlist
   integer, dimension(:), allocatable, intent(out) :: rlist ! row lists
   ! For details of control, info : see derived type descriptions
   type(mc78_control), intent(in) :: control
   integer, intent(out) :: info
   integer, optional, intent(out) :: stat
   integer(long), optional, intent(out) :: nfact ! If present, then on exit
      ! contains the number of entries in L
   integer(long), optional, intent(out) :: nflops ! If present, then on exit
      ! contains the number of floating point operations in factorize.
   integer, dimension(n), optional, intent(inout) :: piv_size ! If
      ! present, then matches matrix order and specifies block pivots. 
      ! piv_size(i) is the number of entries pivots in the block pivot
      ! containing column i.

   integer, dimension(:), allocatable :: cc ! number of entries in each column
   integer :: i
   integer, dimension(:), allocatable :: invp ! inverse permutation of perm
   integer :: j
   integer :: nsvar ! number of supervariables
   integer, dimension(:), allocatable :: parent ! parent of each node in etree
   integer, dimension(:), allocatable :: perm2 ! temporary permutation vector
   integer(pkg_type), dimension(:), allocatable :: ptr ! column pointers for
      ! equivilent matrix
   integer :: realn ! Set to actual number of variables present
   integer, dimension(:), allocatable :: row ! row indices for equivilent matrix
   integer :: st ! stat argument in allocate calls
   integer, dimension(:), allocatable :: svara ! array for supervariables
   integer, dimension(:), allocatable :: sinvp ! inverse permutation of svars
   integer, dimension(:), allocatable :: sperm ! permutation vector of svars
   integer(long) :: sz ! temporary var for size of arrays at allocation

   integer, dimension(:), allocatable :: scc

   ! initialise
   info = 0

   ! Ensure allocatable output arguments are deallocated
   deallocate(sptr, stat=st)
   deallocate(sparent, stat=st)
   deallocate(rptr, stat=st)
   deallocate(rlist, stat=st)

   ! Initialise inverse permutation and check for duplicates
   allocate(invp(n), stat=st)
   if(st.ne.0) goto 490
   do i = 1, n
      j = perm(i)
      invp(j) = i
   end do

   if(present(piv_size)) &
      call convert_to_blk_piv(n, invp, piv_size)

   ! Determine supernodes, build equivilant lwr matrix and find elimination tree
   sz = starts(nelt+1)-1
   if(present(piv_size)) sz = sz + n
   allocate(ptr(n+2), row(sz), svara(n+1), parent(n), stat=st)
   if(st.ne.0) goto 490
   realn = n
   call mc78_elt_equiv_etree(realn, nelt, starts, vars, perm, invp, nsvar, &
      svara, ptr, row, eparent, parent, st, block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Set up permutations of supervariables (initially the idenity)
   allocate(sperm(nsvar), sinvp(nsvar), stat=st)
   if(st.ne.0) goto 490
   sperm(1:nsvar) = (/ (i, i=1,nsvar) /)
   sinvp(1:nsvar) = (/ (i, i=1,nsvar) /)

   ! Postorder tree (modifies perm!)
   call mc78_postorder(nsvar, sperm, sinvp, parent, st, &
      block_pivots=piv_size)
   if(st.ne.0) goto 490

   ! Determine column counts
   allocate(cc(nsvar+1), stat=st)
   if(st.ne.0) goto 490
   call mc78_col_counts(nsvar, ptr, row, sperm, sinvp, parent, cc, st, wt=svara)
   if(st.ne.0) goto 490

   ! Identify supernodes
   allocate(perm2(nsvar), sptr(nsvar+1), sparent(nsvar), scc(nsvar), stat=st)
   if(st.ne.0) goto 490
   call mc78_supernodes(nsvar, nsvar, parent, cc, perm2, nnodes, sptr, &
      sparent, scc, sinvp, control, info, st, wt=svara, block_pivots=piv_size)
   if(info.eq.MC78_ERROR_ALLOC) goto 490
   if(info.lt.0) return

   ! Apply permutation to obtain final elimination order
   call apply_perm(nsvar, perm2, sperm, sinvp, cc, block_pivots=piv_size)

   ! Determine column patterns - keep%nodes(:)%index
   allocate(rptr(nnodes+1), rlist(sum(scc(1:nnodes))), stat=st)
   if(st.ne.0) goto 490
   call mc78_row_lists(nsvar, svara, n, ptr, row, sperm, sinvp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   if(st.ne.0) goto 490

   ! Unmap from supervariables to real variables
   call svar_unmap(n, nsvar, svara, perm, invp, nnodes, sinvp, &
      sptr, st)
   if(st.ne.0) goto 490

   ! Calculate info%num_factor and info%num_flops
   call mc78_stats(nnodes, sptr, scc, nfact=nfact, nflops=nflops)

   if(control%lopt) then
      ! Reorder elimination variables for better cache locality
      if(control%sort) then
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st, sort=.true.)
      else
         call mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, &
            sparent, rptr, rlist, st)
      endif
      if(st.ne.0) goto 490
   else if(control%sort) then
      call dbl_tr_sort(n, nnodes, rptr, rlist, st)
      if(st.ne.0) goto 490
   endif

   ! Adjust eparent with final permuation. Use svara to contain a mapping
   ! from original variables to supernodes
   do i = 1, nnodes
      do j = sptr(i), sptr(i+1)-1
         svara(invp(j)) = i
      end do
   end do
   svara(n+1) = n+1
   do i = 1, nelt
      eparent(i) = svara(eparent(i))
   end do

   if(present(piv_size)) &
      call convert_from_blk_piv(n, perm, piv_size)

   !print *, "n, realn, nelt = ", n, realn, nelt
   !print *, "starts = ", starts
   !print *, "vars = ", vars
   !print *, "nnodes = ", nnodes
   !print *, "sptr = ", sptr(1:nnodes+1)
   !print *, "sparent = ", sparent(1:nnodes)
   !print *, "eparent = ", eparent(:)
   !print *, "rptr = ", rptr(1:nnodes+1)
   !print *, "rlist = ", rlist(1:rptr(nnodes+1)-1)
   !print *, "invp = ", invp
   !print *, "perm = ", perm

   return

   !!!!!!!!!!!!!!!!!!
   ! Error handlers !
   !!!!!!!!!!!!!!!!!!

   490 continue
   info = MC78_ERROR_ALLOC
   if(present(stat)) stat = st
   return
end subroutine mc78_analyse_elemental_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supervariable routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine find supervariables of A using the algorithm of [1].
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
subroutine mc78_supervars_integer(n, ptr, row, perm, invp, nsvar, svar, st)
   integer, intent(inout) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables
   integer, dimension(n), intent(out) :: svar ! number of vars in each svar
   integer, intent(out) :: st

   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occur
   integer(long) :: i
   integer :: j
   integer :: idx ! current index
   integer :: next_sv ! head of free sv linked list
   integer :: nsv ! new supervariable to move j to
   integer :: piv ! current pivot
   integer :: col ! current column of A
   integer :: sv ! current supervariable
   integer :: svc ! temporary holding supervariable count
   integer, dimension(:), allocatable :: sv_new  ! Maps each supervariable to
      ! a new supervariable with which it is associated.
   integer, dimension(:), allocatable :: sv_seen ! Flags whether svariables have
      ! been seen in the current column. sv_seen(j) is set to col when svar j
      ! has been encountered.
   integer, dimension(:), allocatable :: sv_count ! number of variables in sv.

   allocate(sv_new(n+1), sv_seen(n+1), sv_count(n+1), stat=st)
   if(st.ne.0) return

   svar(:) = 1
   sv_count(1) = n
   sv_seen(1) = 0

   ! Setup linked list of free super variables
   next_sv = 2
   do i = 2, n
      sv_seen(i) = i+1
   end do
   sv_seen(n+1) = -1

   ! Determine supervariables using modified Duff and Reid algorithm
   full_rank = .false.
   do col = 1, n
      if(ptr(col+1).ne.ptr(col)) then
         ! If column is not empty, add implicit diagonal entry
         j = col
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! MUST BE the first time that sv has been seen for this
            ! column, so just leave j in sv, and go to next variable.
            ! (Also there can be no other vars in this block pivot)
         else
            ! There is at least one other variable remaining in sv
            ! MUST BE first occurence of sv in the current row/column,
            ! so define a new supervariable and associate it with sv.
            sv_seen(sv) = col
            sv_new(sv) = next_sv
            nsv = next_sv
            next_sv = sv_seen(next_sv)
            sv_new(nsv) = nsv ! avoids problems with duplicates
            sv_seen(nsv) = col
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = 1
            ! This sv cannot be empty as initial sv_count was > 1
         endif
      endif
      do i = ptr(col), ptr(col+1)-1
         j = row(i)
         sv = svar(j)
         if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
            full_rank = full_rank .or. (sv.eq.1)
            ! If so, and this is first time that sv has been seen for this
            ! column, then we can just leave j in sv, and go to next variable.
            if(sv_seen(sv).lt.col) cycle
            ! Otherwise, we have already defined a new supervariable associated
            ! with sv. Move j to this variable, then retire (now empty) sv.
            nsv = sv_new(sv)
            if(sv.eq.nsv) cycle
            svar(j) = nsv
            sv_count(nsv) = sv_count(nsv) + 1
            ! Old sv is now empty, add it to top of free stack
            sv_seen(sv) = next_sv
            next_sv = sv
         else
            ! There is at least one other variable remaining in sv
            if(sv_seen(sv).lt.col) then
               ! this is the first occurence of sv in the current row/column,
               ! so define a new supervariable and associate it with sv.
               sv_seen(sv) = col
               sv_new(sv) = next_sv
               sv_new(next_sv) = next_sv ! avoids problems with duplicates
               next_sv = sv_seen(next_sv)
               sv_count(sv_new(sv)) = 0
               sv_seen(sv_new(sv)) = col
            endif
            ! Now move j from sv to nsv
            nsv = sv_new(sv)
            svar(j) = nsv
            sv_count(sv) = sv_count(sv) - 1
            sv_count(nsv) = sv_count(nsv) + 1
            ! This sv cannot be empty as sv_count was > 1
         endif
      end do
   end do

   ! Note: block pivots do not mix well with supervariables as any significant
   ! number (unless aligned to s.v.s) will demolish any gain from using them.
   ! Converting vlock pivots to s.v.s results in potentially large amount of
   ! unneeded fillin to left of block pivot.
   !! If block pivots are being used, we force all pivots of a block pivot
   !! to be in either the same supervariable, or in supervariables of size 1
   !if(present(block_pivots)) then
   !   piv = 1
   !   do while(piv.le.n)
   !      ! Check if we need to split pivots
   !      split = .false.
   !      sv = svar(piv)
   !      do i = piv+1, piv+block_pivots(piv)
   !         j = invp(i)
   !         if(svar(j).ne.sv) then
   !            split = .true.
   !            exit
   !         endif
   !      end do
   !      ! Do split if required
   !      if(split) then
   !         j = invp(i)
   !         do i = piv, piv+block_pivots(piv)
   !            sv = svar(j)
   !            if(sv_count(sv).eq.1) cycle ! Already a singleton
   !            ! Otherwise create a new sv and move j to it
   !            nsv = next_sv
   !            next_sv = sv_seen(next_sv)
   !            svar(j) = nsv
   !            sv_count(nsv) = 1
   !            sv_count(sv) = sv_count(sv) - 1
   !         end do
   !      endif
   !      piv = piv + block_pivots(piv) + 1
   !   end do
   !endif

   ! Now modify pivot order such that all variables in each supervariable are
   ! consecutive. Do so by iterating over pivots in elimination order. If a
   ! pivot has not already been listed, then order that pivot followed by
   ! any other pivots in that supervariable.

   ! We will build a new inverse permutation in invp, and then find perm
   ! afterwards. First copy invp to perm:
   perm(:) = invp(:)
   ! Next we iterate over the pivots that have not been ordered already
   ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
   ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has been
   ! ordered.
   idx = 1
   nsvar = 0
   do piv = 1, n
      if(sv_seen(piv).gt.n+1) cycle ! already ordered
      ! Record information for supervariable
      sv = svar(perm(piv))
      if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
      nsvar = nsvar + 1
      svc = sv_count(sv)
      sv_new(nsvar) = svc ! store # vars in s.v. to copy to svar later
      j = piv
      ! Find all variables that are members of sv and order them.
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.sv) exit
         end do
         sv_seen(j) = n+2 ! flag as ordered
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
   end do
   ! Push unused variables to end - these are those vars still in s.v. 1
   if(.not.full_rank) then
      svc = sv_count(1)
      ! Find all variables that are members of sv and order them.
      j = 1
      do while(svc.gt.0)
         do j = j, n
            if(svar(perm(j)).eq.1) exit
         end do
         invp(idx) = perm(j)
         idx = idx + 1
         svc = svc - 1
         j = j + 1
      end do
      n = n - sv_count(1)
   end if
   ! Recover perm as inverse of invp
   do piv = 1, n
      perm(invp(piv)) = piv
   end do
   ! sv_new has been used to store number of variables in each svar, copy into
   ! svar where it is returned.
   svar(1:nsvar) = sv_new(1:nsvar)
end subroutine mc78_supervars_integer

!
! This subroutine takes a set of supervariables and compresses the supplied
! matrix using them.
!
! As we would need a full scan of the matrix to calculate the correct size of
! row2, we instead allow the user to make a guess at a good size and return
! an error if this turns out to be incorrect. An upper bound on the required
! size may be obtained by summing the number of entries in the first column of
! each supervariable.
!
! Error returns:
!   MC78_ERROR_ALLOC      Failed to allocate memory
!   MC78_ERROR_ROW_SMALL  row2 too small
subroutine mc78_compress_by_svar_integer(n, ptr, row, invp, nsvar, svar, ptr2, &
      lrow2, row2, info, st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar ! super variables of A
   integer(pkg_type), dimension(nsvar+1), intent(out) :: ptr2
   integer(pkg_type), intent(in) :: lrow2
   integer, dimension(lrow2), intent(out) :: row2
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: piv, svc, sv, col
   integer(pkg_type) :: j, idx
   integer, dimension(:), allocatable :: flag, sv_map

   info = 0 ! by default completed succefully

   allocate(flag(nsvar), sv_map(n), stat=st)
   if(st.ne.0) then
      info = MC78_ERROR_ALLOC
      return
   endif
   flag(:) = 0

   ! Setup sv_map
   piv = 1
   do svc = 1, nsvar
      do piv = piv, piv + svar(svc) - 1
         sv_map( invp(piv) ) = svc
      end do
   end do

   piv = 1
   idx = 1
   do svc = 1, nsvar
      col = invp(piv)
      ptr2(svc) = idx
      do j = ptr(col), ptr(col+1)-1
         sv = sv_map(row(j))
         if(flag(sv).eq.piv) cycle ! Already dealt with this supervariable
         if(idx.gt.lrow2) then
            ! oops, row2 is too small
            info = MC78_ERROR_ROW_SMALL
            return
         endif
         ! Add row entry for this sv
         row2(idx) = sv
         flag(sv) = piv
         idx = idx + 1
      end do
      piv = piv + svar(svc)
   end do
   ptr2(svc) = idx
end subroutine mc78_compress_by_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Elimination tree routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the elimination tree of a PAP^T where A is a
! sparse symmetric matrix stored in compressed sparse column form with
! entries both above and below the diagonal present in the argument matrix.
! P is a permutation stored in order such that order(i) gives the pivot
! position of column i. i.e. order(3) = 5 means that the fifth pivot is
! A_33.
!
! The elimination tree is returned in the array parent. parent(i) gives the
! parent in the elimination tree of pivot i.
!
! The algorithm used is that of Liu [1].
!
! [1] Liu, J. W. 1986. A compact row storage scheme for Cholesky factors using
!     elimination trees. ACM TOMS 12, 2, 127--148.
!
subroutine mc78_etree_integer(n, ptr, row, perm, invp, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer(pkg_type) :: i ! next index into row
   integer :: j ! current entry in row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer :: rowidx ! current column of A = invp(piv)
   integer :: sv ! current supervariable
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   sv = 1
   piv = 1
   do while(piv.le.n)
      !print *, "row ", piv
      rowidx = invp(piv)
      ! Loop over entries in row in lower triangle of PAP^T
      do i = ptr(rowidx), ptr(rowidx+1)-1
         j = perm(row(i))
         if(j.ge.piv) cycle ! not in lower triangle
         !print *, "  entry ", j
         k = j
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         ! Check if we have already done this pivot
         if(vforest(k).eq.piv) cycle 
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
      piv = piv + 1 ! move on to next pivot
   end do
end subroutine mc78_etree_integer

!
! This subroutine identifies supervariables of A using a modified variant of
! the algorithm of Duff and Reid [1]. A lower triangular equivilant matrix
! is returned that is expressed in terms of these supervariables. The grouping
! of variables into supervaribles is returned through a modified pivot order
! and an array specifying the number of variables in each supervariable in
! elimination order. Finally the vector eparent is also returned. This contains
! the variable (in natural numbering) that corresponds to the least pivot in
! each supervariable.
!
! [1] I.S. Duff and J.K. Reid. "Exploiting zeros on the diagonal in the direct
!     solution of indefinite sparse symmetric linear systems".
!     ACM TOMS 22(2), pp227-257. 1996.
!
! Note: If block pivots are present they have priority over supervariables - 
! members of same block pivot must remain in same supervariable. This is
! enforced by moving them all at once.
subroutine mc78_elt_equiv_etree_integer(n, nelt, starts, vars, perm, invp, &
      nsvar, svar, ptr, row, eparent, parent, st, block_pivots)
   integer, intent(inout) :: n ! dimension of system
   integer, intent(in) :: nelt ! number of elements
   integer(pkg_type), dimension(nelt+1), intent(in) :: starts ! variable
      ! pointers of elements
   integer, dimension(starts(nelt+1)-1), intent(in) :: vars ! variables of
      ! elements
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot position
      ! of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, intent(out) :: nsvar ! number of supervariables found
   integer, dimension(n), intent(out) :: svar ! size of each supervariable
   integer(pkg_type), dimension(n+1), intent(out) :: ptr ! column pointers
      ! for equivilant lower triangular form
   integer, dimension(:), intent(out) :: row ! row indices
      ! for equivilant lower triangular form
   integer, dimension(nelt), intent(out) :: eparent ! parent nodes of each
      ! element - i.e. the least pivot in each element
   integer, dimension(n), intent(out) :: parent ! parent(i) is parent of node
      ! i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! Used for
      ! block pivots, see description in analyse phase.

   integer :: csv ! column supervariable in loop
   integer :: elt ! current element
   logical :: full_rank ! flags if supervariable 1 has ever become empty.
      ! If it has not, then the varaibles in s.v. 1 are those that never
      ! occour
   integer(pkg_type) :: i
   integer(pkg_type) :: idx ! Next insert position
   integer :: j
   integer :: k
   integer :: minpiv ! minimum pivot of current element
   integer, dimension(:), allocatable :: mp_head, mp_next ! mp_head and mp_next
      ! store a linked list of elements for which a given variable is the
      ! minimum pivot.
   integer :: next_sv ! Top of stack of free supervariables (stored as a linked
      ! list in unused part of sv_seen)
   integer :: orign ! original system dimension
   integer :: nsv ! temporary variable storing new supervariable to move var to
   integer :: piv ! current pivot
   integer :: sv ! current supervariable
   integer :: svc ! temporary variable storing supervariable count remaining
   integer, dimension(:), allocatable :: sv_count ! sv_count(s) is the number
      ! of variables in supervariable s.
   integer, dimension(:), allocatable :: sv_map ! sv_map(v) is the current
      ! supervariable to which variable v belongs.
   integer, dimension(:), allocatable :: sv_new ! sv_map(s) is new
      ! supervariable for variables currently in supervariable s.
   integer, dimension(:), allocatable :: sv_seen ! sv_seen(s) is used to flag
      ! if supervariable s has been found in the current element before.
      ! In addition the part corresponding to unused supervariables is used
      ! to store a stack (as a linked list) of empty supervariables.
   integer(pkg_type), dimension(:), allocatable :: uprptr ! column pointers for
      ! upper triangular equivilant form.
   integer, dimension(:), allocatable :: uprrow ! row indices for upper
      ! triangular equivilant form.
   logical :: used ! flag if a variable has been used

   ! Initialise supervariable representation
   allocate(sv_new(max(nelt+1,n+1)), sv_seen(max(nelt,n)+1), &
      sv_map(max(nelt,n)+1), sv_count(max(nelt,n)+1), stat=st)
   if(st.ne.0) return
   sv_map(:) = 1 ! All vars are intially in supervariable 1
   sv_count(1) = n ! ... which thus has all variables
   sv_seen(1) = 0 ! Flag supervariable 1 as unseen on first iteration
   orign = n

   if(present(block_pivots)) then
      ! Do not mix supervariables and block pivots

      ! Still need to determine minimum pivots and rank
      sv_seen(:) = 0
      do elt = 1, nelt
         minpiv = n+1
         do i = starts(elt), starts(elt+1)-1
            j = vars(i)
            ! Mark variable as used
            sv_seen(j) = 1
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) eparent(elt) = invp(minpiv)
      end do

      ! Build invp that pushes unsued vars to the end. Be careful of unused vars
      ! that are in fact part of block pivots, split them out.
      perm(:) = invp(:)
      piv = 1
      j = 1
      ! Handle variables that are actually used
      do while(piv.le.n)
         used = .true.
         do i = piv, n
            used = used .and. (sv_seen(perm(i)).eq.1)
            if(block_pivots(i).ge.2) exit ! end of block pivot
         end do

         if(used) then
            ! Block pivot is entirely composed of used variables
            do piv = piv, i
               invp(j) = perm(piv)
               sv_seen(perm(piv)) = 2
               j = j + 1
            end do
         else
            ! Block pivot has some unused variables in it
            k = 0
            do piv = piv, i
               if(sv_seen(perm(piv)).eq.1) then
                  invp(j) = perm(piv)
                  sv_seen(perm(piv)) = 2
                  j = j + 1
                  if(k.eq.0) then
                     ! This is the new start of the block pivot
                     select case(block_pivots(piv))
                     case(0) ! was in the middle. now a start
                        block_pivots(piv) = 1
                     case(2) ! was the end. now a 1x1
                        block_pivots(piv) = 3
                     end select
                  endif
                  k = piv
               endif
            end do
            if(k.ne.0) then
               ! The was at least one used variable in the block pivot
               select case(block_pivots(k))
               case(0) ! was the middle. now an end
                  block_pivots(k) = 2
               case(1) ! was the start. now a 1x1
                  block_pivots(k) = 3
               end select
            endif
         endif
         piv = i + 1
      end do
      ! Handle unused variables; build supervariables
      nsvar = 0
      do piv = 1, n
         i = perm(piv)
         if(sv_seen(i).eq.0) then
            invp(j) = i
            j = j + 1
            block_pivots(piv) = 3 ! Force to 1x1
         else
            nsvar = nsvar + 1
            svar(nsvar) = 1
            sv_new(i) = nsvar
         endif
      end do
      ! Map block_pivots in original variable order into sv_map
      do i = 1, n
         sv_map(perm(i)) = block_pivots(i)
      end do
      ! Map sv_map in new pivot order back into block_pivots
      do i = 1, n
         block_pivots(i) = sv_map(invp(i))
      end do
      ! Reestablish perm
      do i = 1, n
         perm(invp(i)) = i
      end do
   else
      ! Setup linked list of free supervariables - we utilise the unused part
      ! of sv_seen for this
      next_sv = 2
      do i = 2, n
         sv_seen(i) = i+1
      end do
      sv_seen(n+1) = -1
      
      ! Determine supervariables. At the same time find the least pivot
      ! associated with each element.
      nsvar = 1
      full_rank = .false.
      do elt = 1, nelt
         minpiv = n+1
         do i = starts(elt), starts(elt+1)-1
            j = vars(i)
            ! Determine minimum pivot
            minpiv = min(minpiv, perm(j))
            sv = sv_map(j)
            if(sv_count(sv).eq.1) then ! Are we only (remaining) var in sv
               full_rank = full_rank .or. (sv.eq.1)
               ! If so, and this is first time that sv has been seen for this
               ! element, then we can just leave j in sv, and go to next
               ! variable.
               if(sv_seen(sv).lt.elt) cycle
               ! Otherwise, we have already defined a new supervariable
               ! associated with sv. Move j to this variable, then retire (now
               ! empty) sv.
               ! Note: as only var in sv, cannot have fellows in block pivot
               nsv = sv_new(sv)
               if(sv.eq.nsv) cycle  ! don't delete a variable because of a
                                    ! duplicate
               sv_map(j) = nsv
               sv_count(nsv) = sv_count(nsv) + 1
               ! Old sv is now empty, add it to top of free stack
               sv_seen(sv) = next_sv
               next_sv = sv
               nsvar = nsvar - 1
            else
               ! There is at least one other variable remaining in sv
               if(sv_seen(sv).lt.elt) then
                  ! this is the first occurence of sv in the current element,
                  ! so define a new supervariable and associate it with sv.
                  sv_seen(sv) = elt
                  sv_new(sv) = next_sv
                  sv_new(next_sv) = next_sv  ! ensure we are tolerant of
                                             ! duplicates
                  next_sv = sv_seen(next_sv)
                  sv_count(sv_new(sv)) = 0
                  sv_seen(sv_new(sv)) = elt
                  nsvar = nsvar + 1
               endif
               ! Now move j from sv to nsv
               nsv = sv_new(sv)
               sv_map(j) = nsv
               sv_count(sv) = sv_count(sv) - 1
               sv_count(nsv) = sv_count(nsv) + 1
               ! We know sv can't be empty, so it doesn't need adding to free
               ! stack
            endif
         end do
         ! Store the minimum pivot as original variable (avoids messy remapping)
         if(minpiv.le.n) then
            eparent(elt) = invp(minpiv)
         else
            eparent(elt) = n+1
         endif
      end do

      ! Now modify pivot order such that all variables in each supervariable are
      ! consecutive. Do so by iterating over pivots in elimination order. If a
      ! pivot has not already been listed, then order that pivot followed by
      ! any other pivots in that supervariable.

      ! We will build a new inverse permutation in invp, and then find perm
      ! afterwards. First copy invp to perm:
      perm(:) = invp(:)
      ! Next we iterate over the pivots that have not been ordered already
      ! Note: as we begin, all entries of sv_seen are less than or equal to n+1
      ! hence we can use <=n+1 or >n+1 as a flag to indicate if a variable has
      ! been ordered.
      idx = 1
      nsvar = 0
      do piv = 1, n
         if(sv_seen(piv).gt.n+1) cycle ! already ordered
         ! Record information for supervariable
         sv = sv_map(perm(piv))
         if(.not.full_rank .and. sv.eq.1) cycle ! Don't touch unused vars
         svc = sv_count(sv)
         nsvar = nsvar + 1
         svar(nsvar) = svc
         ! Find all variables that are members of sv and order them.
         j = piv
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.sv) exit
            end do
            sv_seen(j) = n+2 ! flag as ordered
            sv_new(perm(j)) = nsvar ! new mapping to sv
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
      end do
      sv_new(n+1) = nsvar+1
      ! Push unused variables to end - these are those vars still in s.v. 1
      if(.not.full_rank) then
         svc = sv_count(1)
         ! Find all variables that are members of sv and order them.
         j = 1
         do while(svc.gt.0)
            do j = j, n
               if(sv_map(perm(j)).eq.1) exit
            end do
            invp(idx) = perm(j)
            idx = idx + 1
            svc = svc - 1
            j = j + 1
         end do
         n = n - sv_count(1)
      end if
      ! Recover perm as inverse of invp
      do piv = 1, n
         perm(invp(piv)) = piv
      end do
   endif

   ! build linked lists by supervariable
   allocate(mp_head(nsvar+1), mp_next(nelt), stat=st)
   if(st.ne.0) return
   mp_head(:) = -1
   do elt = 1, nelt
      if(eparent(elt).gt.orign) cycle
      minpiv = sv_new(eparent(elt))
      if(present(block_pivots) .and. minpiv.ne.1) then
         do while(block_pivots(minpiv-1).lt.2)
            minpiv = minpiv - 1
            if(minpiv.eq.1) exit
         end do
      endif
      ! Store element in linked list for minpiv
      mp_next(elt) = mp_head(minpiv)
      mp_head(minpiv) = elt
   end do

   ! Iterate over columns in pivot order, storing the lower triangular
   ! equivilant matrix as we go. At the same time, build the column counts for
   ! the upper triangle in uprptr, but offset by 2 (ie uprptr(i+2) for col i).
   ! Observe that all the pivots associated with the supervariable
   ! to which minpiv belongs _must_ appear in each element that minpiv does.
   ! Note: This only generates the lower triangular part of the matrix!
   allocate(uprptr(nsvar+2), stat=st)
   if(st.ne.0) return
   uprptr(:) = 0
   sv_seen(:) = 0
   idx = 1
   ptr(:) = -1
   do csv = 1, nsvar
      elt = mp_head(csv)
      ptr(csv) = idx
      sv_seen(csv) = csv ! Mark diagonal as seen, as it is implicit.
      do while(elt.ne.-1)
         do i = starts(elt), starts(elt+1)-1
            sv = sv_new(vars(i))
            ! Skip this sv if it is already included (or is implicit)
            if(sv_seen(sv).ge.csv) cycle
            sv_seen(sv) = csv ! Mark as seen
            ! If we can't skip it, then add entry (sv,csv) to lwr matrix
            row(idx) = sv
            idx = idx + 1
            ! Add count in upper triangle for (csv,sv)
            uprptr(sv+2) = uprptr(sv+2) + 1
         end do
         ! Move on to next element for which this is the minimum pivot
         elt = mp_next(elt)
      end do
      if(present(block_pivots) .and. csv.ne.nsvar) then
         ! Add entry (csv+1,csv) to ensure elimination tree correct
         if(block_pivots(csv).lt.2 .and. sv_seen(csv+1).ne.csv) then
            sv_seen(csv+1) = csv
            row(idx) = csv+1
            idx = idx + 1
            ! Add count in upper triangle for (csv, csv+1)
            uprptr(csv+1+2) = uprptr(csv+1+2) + 1
         endif
      endif
   end do
   ptr(nsvar+1) = idx

   ! Build upper form - work out column start for col i in uprptr(i+1)
   uprptr(1:2) = 1
   do i = 1, nsvar
      uprptr(i+2) = uprptr(i+1) + uprptr(i+2)
   end do

   ! Now iterate over lwr form, droppping entries into upr form
   allocate(uprrow(uprptr(nsvar+2)), stat=st)
   if(st.ne.0) return
   do csv = 1, nsvar
      do i = ptr(csv), ptr(csv+1) - 1
         sv = row(i)
         uprrow(uprptr(sv+1)) = csv
         uprptr(sv+1) = uprptr(sv+1) + 1
      end do
   end do

   ! Now determine supervariable elimination tree
   call etree_no_perm(nsvar, uprptr, uprrow, parent, st)
   if(st.ne.0) return
end subroutine mc78_elt_equiv_etree_integer

!
! Specialised version of mc78_etree that assumes elimination order is identity
!
subroutine etree_no_perm(n, ptr, row, parent, st)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(out) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls

   integer :: i ! next index into row
   integer :: k ! current ancestor
   integer :: l ! next ancestor
   integer :: piv ! current pivot
   integer, dimension(:), allocatable :: vforest ! virtual forest, used for
      ! path compression (shortcuts to top of each tree)

   ! Allocate virtual forest and initialise it
   allocate(vforest(n), stat=st)
   if(st.ne.0) return
   vforest(:) = n+1

   ! Loop over rows of A in pivot order
   do piv = 1, n
      ! Loop over entries in row in lower triangle of PAP^T
      do i = ptr(piv), ptr(piv+1)-1
         k = row(i)
         do while(vforest(k).lt.piv)
            l = vforest(k)
            vforest(k) = piv
            k = l
         end do
         if(vforest(k).eq.piv) cycle ! Already done from here, don't overwrite
         parent(k) = piv
         vforest(k) = piv
      end do
      parent(piv) = n + 1 ! set to be a root if not overwritten
   end do
end subroutine etree_no_perm

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_std(n, perm, invp, parent, st, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_std

!
! This subroutine will postorder the elimination tree. That is to say it will
! reorder the nodes of the tree such that they are in depth-first search order.
!
! This is done by performing a depth-first search to identify mapping from the
! original pivot order to the new one. This map is then applied to order, invp
! and parent to enact the relabelling.
!
subroutine mc78_postorder_detect(n, realn, ptr, perm, invp, parent, st, &
      block_pivots)
   integer, intent(in) :: n
   integer, intent(out) :: realn
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(n), intent(inout) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(inout) :: invp ! inverse of perm
   integer, dimension(n), intent(inout) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(n), optional, intent(inout) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer :: i
   integer :: id
   integer :: j
   integer, dimension(:), allocatable :: map ! mapping from original pivot
      ! order to new one
   integer :: node
   integer :: shead ! pointer to top of stack
   integer, dimension(:), allocatable :: stack ! stack for depth first search

   realn = n

   !
   ! Build linked lists of children for each node
   !
   allocate(chead(n+1), cnext(n+1), stat=st)
   if(st.ne.0) return
   chead(:) = -1 ! no parent if necessary
   do i = n, 1, -1 ! do in reverse order so they come off in original order
      j = parent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search to build map
   !
   allocate(map(n+1), stack(n), stat=st)
   if(st.ne.0) return
   ! Place virtual root on top of stack
   shead = 1
   stack(shead) = n+1
   id = n + 1 ! next node id
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      map(node) = id
      id = id - 1

      ! Place all its children on the stack such that the last child is
      ! at the top of the stack and first child closest to the bottom
      if(node.eq.n+1) then
         ! Virtual root node, detect children with no entries at same time
         ! placing those that are empty at the top of the stack
         ! First do those which are proper roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).eq.0) then
               i = cnext(i)
               cycle
            endif
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
         ! Second do those which are null roots
         i = chead(node)
         do while(i.ne.-1)
            if(ptr(invp(i)+1)-ptr(invp(i)).ne.0) then
               i = cnext(i)
               cycle
            endif
            realn = realn - 1
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      else ! A normal node
         i = chead(node)
         do while(i.ne.-1)
            shead = shead + 1
            stack(shead) = i
            i = cnext(i)
         end do
      endif
   end do

   !
   ! Apply map to perm, invp and parent (and block_pivots if present)
   !

   ! invp is straight forward, use stack as a temporary
   stack(1:n) = invp(1:n)
   do i = 1, n
      j = map(i)
      invp(j) = stack(i)
   end do

   ! perm can be easily done as the inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! parent is done in two stages. The first copies it to stack and permutes
   ! parent(i), but not the locations. i.e. if 1 is a parent of 3, and
   ! map(1)=2 and map(3)=4, then the first stage sets stack(1) = 4.
   ! The second stage then permutes the entries of map back into parent
   do i = 1, n
      stack(i) = map(parent(i))
   end do
   do i = 1, n
      parent(map(i)) = stack(i)
   end do

   ! permute block_pivots if required
   if(present(block_pivots)) then
      stack(1:n) = block_pivots(1:n)
      do i = 1, n
         block_pivots(map(i)) = stack(i)
      end do
   endif
end subroutine mc78_postorder_detect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Column count routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines column counts given the elimination tree and
! pattern of the matrix PAP^T.
!
! The algorithm is a specialisation of that given by Gilbert, Ng and Peyton [1],
! to only determine column counts. It is also described in Section 4.4 "Row
! counts" of [2].
!
! The essential technique is to determine the net number of entries introduced
! at a node (the "weight" in [1]). This is composed over the following terms:
!  wt[i] = [ - #children of node i
!            - #common indices between children
!            + #additional "new" row indices from column of A ]
!
! The clever part of this algorithm is how to determine the number of common
! indices between the children. This is accomplished by storing the last column
! at which an index was encountered, and a partial elimination tree. This
! partial elimination tree consists of all nodes processed so far, plus their
! parents. As we have a postorder on the tree, the current top of the tree
! containing node i is the least common ancestor of node i and the current node.
! We then observe that the first time an index will be double counted is at the
! least common ancestor of the current node and the last node where it was
! encountered.
!
! [1] Gilbert, Ng, Peyton, "An efficient algorithm to compute row and column
!     counts for sparse Cholesky factorization", SIMAX 15(4) 1994.
!
! [2] Tim Davis's book "Direct Methods for Sparse Linear Systems", SIAM 2006.
!
subroutine mc78_col_counts_integer(n, ptr, row, perm, invp, parent, cc, st, wt)
   integer, intent(in) :: n ! dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! column pointers of A
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices of A
   integer, dimension(n), intent(in) :: perm ! perm(i) is the pivot
      ! position of column i
   integer, dimension(n), intent(in) :: invp ! inverse of perm
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of pivot i in the elimination tree
   integer, dimension(n+1), intent(out) :: cc ! On exit, cc(i) is the
      ! number of entries in the lower triangular part of L (includes diagonal)
      ! for the column containing pivot i. For most of the routine however, it
      ! is used as a work space to track the net number of entries appearing
      ! for the first time at node i of the elimination tree (this may be
      ! negative).
   integer, intent(out) :: st ! stat parmeter for allocate calls
   integer, dimension(:), optional, intent(in) :: wt ! weights (eg number of
      ! variables in each supervariable)
   
   integer :: col ! column of matrix associated with piv
   integer, dimension(:), allocatable :: first ! first descendants
   integer(pkg_type) :: i ! loop variable
   integer :: totalwt
   integer, dimension(:), allocatable :: last_nbr ! previous neighbour
   integer, dimension(:), allocatable :: last_p ! previous p?
   integer :: par ! parent node of piv
   integer :: piv ! current pivot
   integer :: pp ! last pivot where u was encountered
   integer :: lca ! least common ancestor of piv and pp
   integer :: u ! current entry in column col
   integer :: uwt ! weight of u
   integer, dimension(:), allocatable :: vforest ! virtual forest

   !
   ! Determine first descendants, and set cc = 1 for leaves and cc = 0 for
   ! non-leaves.
   !
   allocate(first(n+1), stat=st)
   if(st.ne.0) return
   do i = 1, n+1
      first(i) = i
   end do
   if(present(wt)) then
      totalwt = 0 ! Find sum of weights so we can determine non-physical value
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = wt(invp(i))
         else
            cc(i) = 0
         endif
         totalwt = totalwt + wt(invp(i))
      end do
      cc(n+1) = totalwt + 1 ! Set to non-physical value
   else
      do i = 1, n
         par = parent(i)
         first(par) = min(first(i), first(par)) ! first descendant
         if(first(i).eq.i) then ! is it a leaf or not?
            cc(i) = 1
         else
            cc(i) = 0
         endif
      end do
      cc(n+1) = n + 1 ! Set to non-physical value
   endif

   !
   ! We store the partial elimination trees in a virtual forest. It is
   ! initialised such that each node is in its own tree to begin with.
   !
   allocate(vforest(n+1), stat=st)
   if(st.ne.0) return
   vforest(:) = 0

   !
   ! Initialise previous pivot and neightbour arrays to indicate no previous
   ! pivot or neightbour.
   !
   allocate(last_p(n+1), last_nbr(n+1), stat=st)
   if(st.ne.0) return
   last_p(:) = 0
   last_nbr(:) = 0

   !
   ! Determine cc(i), the number of net new entries to pass up tree from
   ! node i.
   !
   do piv = 1, n
      ! Loop over entries in column below the diagonal
      col = invp(piv)
      do i = ptr(col), ptr(col+1)-1
         u = perm(row(i))
         if(u.le.piv) cycle ! not in lower triangular part

         ! Check if entry has been seen by a descendant of this pivot, if
         ! so we skip the tests that would first add one to the current
         ! pivot's weight before then subtracting it again.
         if(first(piv).gt.last_nbr(u)) then
            ! Count new entry in current column
            uwt = 1
            if(present(wt)) uwt = wt(invp(u))
            cc(piv) = cc(piv) + uwt

            ! Determine least common ancestor of piv and the node at which
            ! u was last encountred
            pp = last_p(u)
            if(pp.ne.0) then
               ! u has been seen before, find top of partial elimination
               ! tree for node pp
               lca = FIND(vforest, pp)
               ! prevent double counting of u at node lca
               cc(lca) = cc(lca) - uwt
            endif

            ! Update last as u has now been seen at piv.
            last_p(u) = piv
         endif

         ! Record last neighbour of u so we can determine if it has been
         ! seen in this subtree before
         last_nbr(u) = piv
      end do
      ! Pass uneliminated variables up to parent
      par = parent(piv)
      if(present(wt)) then
         cc(par) = cc(par) + cc(piv) - wt(invp(piv))
      else
         cc(par) = cc(par) + cc(piv) - 1
      endif

      ! place the parent of piv into the same partial elimination tree as piv
      vforest(piv) = par ! operation "UNION" from [1]
   end do
end subroutine mc78_col_counts_integer

! Return top most element of tree containing u.
! Implements path compression to speed up subsequent searches.
integer function FIND(vforest, u)
   integer, dimension(:), intent(inout) :: vforest
   integer, intent(in) :: u

   integer :: current, prev

   prev = -1
   current = u
   do while(vforest(current).ne.0)
      prev = current
      current = vforest(current)
      if(vforest(current).ne.0) vforest(prev) = vforest(current)
   end do

   FIND = current
end function FIND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Supernode amalgamation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine identifies (relaxed) supernodes from the elimination tree
! and column counts.
!
! A node, u, and its parent, v, are merged if:
! (a) No new fill-in is introduced i.e. cc(v) = cc(u)-1
! (b) The number of columns in both u and v is less than nemin
!
! Note: assembly tree must be POSTORDERED on output
subroutine mc78_supernodes(n, realn, parent, cc, sperm, nnodes, sptr, sparent, &
      scc, invp, control, info, st, wt, block_pivots)
   integer, intent(in) :: n
   integer, intent(in) :: realn
   integer, dimension(n), intent(in) :: parent ! parent(i) is the
      ! parent of supernode i in the elimination/assembly tree. 
   integer, dimension(n), intent(in) :: cc ! cc(i) is the column count
      ! of supernode i, including elements eliminated at supernode i.
   integer, dimension(n), intent(out) :: sperm ! on exit contains a permutation
      ! from pivot order to a new pivot order with contigous supernodes
   integer, intent(out) :: nnodes ! number of supernodes
   integer, dimension(n+1), intent(out) :: sptr
   integer, dimension(n), intent(out) :: sparent
   integer, dimension(n), intent(out) :: scc
   integer, dimension(n), intent(in) :: invp
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st ! stat paremter from allocate calls
   integer, dimension(n), optional, intent(in) :: wt ! weights (number of vars
      ! in each initial node)
   integer, dimension(n), optional, intent(in) :: block_pivots ! If
      ! present, then matches pivot order and specifies block pivots. 
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot

   integer :: i, j, k
   integer :: flag
   integer, dimension(:), allocatable :: height ! used to track height of tree
   logical, dimension(:), allocatable :: mark ! flag array for nodes to finalise
   integer, dimension(:), allocatable :: map ! map vertex idx -> supernode idx
   integer, dimension(:), allocatable :: nelim ! number of eliminated variables
   integer, dimension(:), allocatable :: nvert ! number of elimd supervariables
   integer :: node
   integer, dimension(:), allocatable :: npar ! temporary array of snode pars
   integer :: par ! parent of current node
   integer :: shead ! current head of stack
   integer, dimension(:), allocatable :: stack ! used to navigate tree
   integer :: v
   integer, dimension(:), allocatable :: vhead ! heads of vertex linked lists
   integer, dimension(:), allocatable :: vnext ! next element in linked lists
   integer(long), dimension(:), allocatable :: ezero ! number of explicit zeros
   integer, dimension(:), allocatable :: chead ! chead(i) is first child of i
   integer, dimension(:), allocatable :: cnext ! cnext(i) is next child of i
   integer, dimension(:), allocatable :: child
   integer :: nchild
   integer :: start ! First pivot in block pivot
   integer :: totalwt ! sum of weights

   !
   ! Initialise supernode representation
   !
   allocate(nelim(n+1), nvert(n+1), vhead(n+1), vnext(n+1), stack(n), &
      height(n+1), mark(n), stat=st)
   if(st.ne.0) goto 490
   vnext(:) = -1
   vhead(:) = -1
   height(:) = 1

   ! Initialise number of variables in each node
   if(present(wt)) then
      totalwt = 0
      do i = 1, n
         nelim(i) = wt(invp(i))
         totalwt = totalwt + wt(invp(i))
      end do
   else ! All nodes initially contain a single variable
      nelim(:) = 1
      totalwt = n
   endif
   nvert(:) = 1

   allocate(map(n+1), npar(n+1), ezero(n+1), stat=st)
   if(st.ne.0) goto 490

   ezero(:) = 0 ! Initially no explicit zeros
   ezero(n+1) = huge(ezero) ! ensure root is not merged

   ! Ensure virtual root never gets amlgamated
   nelim(n+1) = totalwt+1 + control%nemin

   !
   ! Build child linked lists for nodes; merge block pivots if needed
   !
   allocate(chead(n+1), cnext(n+1), child(n), stat=st)
   if(st.ne.0) goto 490
   if(present(block_pivots)) then
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         if(block_pivots(i).lt.2) cycle
         j = parent(i)
         if(j.ne.n+1) then
            do while(block_pivots(j).lt.2)
               j = parent(j)
            end do
         end if
         cnext(i) = chead(j)
         chead(j) = i
      end do
   else
      chead(:) = -1 ! no parent if necessary
      do i = realn, 1, -1 ! do in reverse order so come off in original order
         j = parent(i)
         cnext(i) = chead(j)
         chead(j) = i
      end do
   endif

   !
   ! Merge supernodes.
   !
   v = 1
   nnodes = 0
   start=n+2
   do par = 1, n+1
      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).lt.2) then
            if(start.ge.n+1) start = par
            cycle
         endif
         ! Merge pivots start to par, but don't add to vertex list (yet)

         do node = start, par-1
            ! Add together eliminated variables
            nelim(par) = nelim(par) + nelim(node)
            nvert(par) = nvert(par) + nvert(node)

            ! nodes have same height
            height(par) = max(height(par), height(node))
         end do
      endif

      nchild = 0
      node = chead(par)
      do while(node.ne.-1)
         nchild = nchild + 1
         child(nchild) = node
         node = cnext(node)
      end do
      call sort_by_val(nchild, child, cc, st)
      if(st.ne.0) goto 490

      do j = 1, nchild
         node = child(j)
         if(do_merge(node, par, nelim, cc, ezero, control, invp, flag, wt)) then
            ! Merge contents of node into par. Delete node.
            call merge_nodes(node, par, nelim, nvert, vhead, vnext, height, &
               ezero, cc)
            mark(node) = .false.
         else
            mark(node) = .true.
         endif
      end do

      if(present(block_pivots) .and. par.lt.n+1) then
         if(block_pivots(par).ge.2) then
            ! Add vertices start to par-1 into par
            do node = start, par-1
               vnext(node) = vhead(par)
               vhead(par) = node
               mark(node) = .false.
            end do
            start = n+2
         endif
      endif
   end do

   if(flag.ne.0) then
      if(control%unit_error.gt.0) write(control%unit_error, "(a)") &
         "MC78 Internal Error: Unrecognised amalgamation heuristic."
      info = MC78_ERROR_UNKNOWN
      return
   endif

   do node = 1, realn
      if(.not.mark(node)) cycle
      ! Node not merged, now a complete supernode

      ! Record start of supernode
      nnodes = nnodes + 1
      sptr(nnodes) = v
      npar(nnodes) = parent(node)
      if(present(wt)) then
         scc(nnodes) = cc(node) + nelim(node) - wt(invp(node))
      else
         scc(nnodes) = cc(node) + nelim(node) - 1
      endif

      ! Record height in tree of parent vertices
      height(parent(node)) = max(height(parent(node)), height(node) + 1)

      ! Determine last vertex of node so we can number backwards
      v = v + nvert(node)
      k = v

      ! Loop over member vertices of node and number them
      shead = 1
      stack(shead) = node
      do while(shead.gt.0)
         i = stack(shead)
         shead = shead - 1

         ! Order current vertex
         k = k - 1
         sperm(i) = k
         map(i) = nnodes

         ! Stack successor, if any
         if(vnext(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vnext(i)
         endif

         ! Descend into tree rooted at i
         if(vhead(i).ne.-1) then
            shead = shead + 1
            stack(shead) = vhead(i)
         endif
      end do
   end do
   sptr(nnodes+1) = v ! Record end of final supernode
   map(n+1) = nnodes + 1 ! virtual root vertex maps to virtual root sn
   npar(nnodes+1) = n + 1

   ! Handle permutation of empty columns
   do i = realn+1, n
      sperm(i) = i
   end do

   ! Allocate arrays for return and copy data into them correctly
   do node = 1, nnodes
      par = npar(node) ! parent /vertex/ of supernode
      par = map(par)   ! parent /node/   of supernode
      sparent(node) = par ! store parent
   end do

   return

   490 continue
   info = MC78_ERROR_ALLOC
   return
end subroutine mc78_supernodes

!
! Sort n items labelled by idx into decreasing order of val(idx(i))
!
recursive subroutine sort_by_val(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: ice_idx, ice_val, ik_idx, ik_val
   integer :: klo,kor,k,kdummy

   st = 0

   if(n.ge.minsz_ms) then
      call sort_by_val_ms(n, idx, val, st)
   else
      klo = 2
      kor = n
      do kdummy = klo,n
         ! items kor, kor+1, .... ,nchild are in order
         ice_idx = idx(kor-1)
         ice_val = val(ice_idx)
         do k = kor,n
            ik_idx = idx(k)
            ik_val = val(ik_idx)
            if (ice_val >= ik_val) exit
            idx(k-1) = ik_idx
         end do
         idx(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine sort_by_val

! Sort n items labelled by idx into decreasing order of val(idx(i))
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine sort_by_val_ms(n, idx, val, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: idx
   integer, dimension(:), intent(in) :: val
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call sort_by_val(n, idx, val, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call sort_by_val_ms(mid, idx(1:mid), val, st)
   if(st.ne.0) return
   call sort_by_val_ms(n - mid, idx(mid+1:n), val, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = idx(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   jj2 = val(jj)
   kk = idx(k)
   kk2 = val(kk)
   do i = 1, n
      if(jj2.ge.kk2) then
         idx(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         jj2 = val(jj)
      else
         idx(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = idx(k)
         kk2 = val(kk)
      endif
   end do
   if(j.le.mid) idx(i+1:n) = work(j:mid)
end subroutine sort_by_val_ms

!
! Return .true. if we should merge node and par, .false. if we should not
!
logical function do_merge(node, par, nelim, cc, ezero, control, invp, info, wt)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(in) :: nelim
   integer, dimension(:), intent(in) :: cc
   integer(long), dimension(:), intent(in) :: ezero
   type(mc78_control), intent(in) :: control
   integer, dimension(:), intent(in) :: invp
   integer, intent(out) :: info
   integer, dimension(:), optional, intent(in) :: wt

   real(dp) :: z, ne

   info = 0

   if(ezero(par).eq.huge(ezero)) then
      do_merge = .false.
      return
   endif

   select case(control%heuristic)
   case(1)
      !
      ! HSL_MA77 style nemin
      !
      if(present(wt)) then
         do_merge = (cc(par).eq.cc(node)-wt(invp(node)) .and. &
            nelim(par).eq.wt(invp(par))) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      else
         do_merge = (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
            (nelim(par).lt.control%nemin .and. nelim(node).lt.control%nemin)
      endif
   case(2)
      !
      ! CHOLMOD style nrelax/zrelax
      !
      ! FIXME: currently assumes nodes are square, not trapezoidal

      ! calculate number of non-zeros in new node
      z = ezero(par) + ezero(node) + &
         (cc(par)-1+nelim(par) - cc(node)+1) * nelim(par)
      ! find this as a fraction of total non-zeros in new node
      ne = nelim(par) + nelim(node)
      z = z / ( (cc(par)-1+ne)*ne )

      do_merge = (ne .le. control%nrelax(1)) .or. &
         (cc(par).eq.cc(node)-1 .and. nelim(par).eq.1) .or. &
         (ne .le. control%nrelax(2) .and. z .lt. control%zrelax(1)) .or. &
         (ne .le. control%nrelax(3) .and. z .lt. control%zrelax(2)) .or. &
         (z .lt. control%zrelax(3))
   case default
      ! Note: This bit of code should NEVER execute
      do_merge = .false.
      info = MC78_ERROR_UNKNOWN
   end select
end function do_merge

!
! This subroutine merges node with its parent, deleting node in the process.
!
subroutine merge_nodes(node, par, nelim, nvert, vhead, vnext, height, ezero, cc)
   integer, intent(in) :: node ! node to merge and delete
   integer, intent(in) :: par ! parent to merge into
   integer, dimension(:), intent(inout) :: nelim
   integer, dimension(:), intent(inout) :: nvert
   integer, dimension(:), intent(inout) :: vhead
   integer, dimension(:), intent(inout) :: vnext
   integer, dimension(:), intent(inout) :: height
   integer(long), dimension(:), intent(inout) :: ezero
   integer, dimension(:), intent(in) :: cc

   ! Add node to list of children merged into par
   vnext(node) = vhead(par)
   vhead(par) = node

   ! Work out number of explicit zeros in new node
   ! FIXME: probably wrong now with weights and block pivots
   ezero(par) = ezero(par) + ezero(node) + &
      (cc(par)-1+nelim(par) - cc(node) + 1_long) * nelim(par)

   ! Add together eliminated variables
   nelim(par) = nelim(par) + nelim(node)
   nvert(par) = nvert(par) + nvert(node)

   ! nodes have same height
   height(par) = max(height(par), height(node))
end subroutine merge_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Statistics routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine merely calculates interesting statistics
!
subroutine mc78_stats(nnodes, sptr, scc, nfact, nflops)
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), optional, intent(out) :: nfact
   integer(long), optional, intent(out) :: nflops

   integer :: j
   integer :: m ! number of entries in retangular part of ndoe
   integer :: nelim ! width of node
   integer :: node ! current node of assembly tree
   integer(long) :: r_nfact, r_nflops

   if(.not.present(nfact) .and. .not.present(nflops)) return ! nothing to do

   r_nfact = 0
   r_nflops = 0
   do node = 1, nnodes
      nelim = sptr(node+1) - sptr(node)
      m = scc(node) - nelim

      ! number of entries
      r_nfact = r_nfact + (nelim * (nelim+1)) / 2 ! triangular block
      r_nfact = r_nfact + nelim * m ! below triangular block

      ! flops
      do j = 1, nelim
         r_nflops = r_nflops + (m+j)**2
      end do
   end do

   if(present(nfact)) nfact = r_nfact
   if(present(nflops)) nflops = r_nflops

   !print *, "n = ", n
   !print *, "nnodes = ", nnodes
   !print *, "nfact = ", nfact
   !print *, "sum cc=", sum(cc(1:n))
   !print *, "nflops = ", nflops
end subroutine mc78_stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Row list routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This subroutine determines the row indices for each supernode
!
subroutine mc78_row_lists_nosvar_integer(n, ptr, row, perm, invp, nnodes, &
      sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: n
   integer(pkg_type), dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: nelim ! number of variables eliminated at node
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child
   integer :: sz ! number of rows in node

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), stat=st)
   if(st.ne.0) return
   seen(:) = 0
   chead(:) = -1

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      nelim = sptr(node+1) - sptr(node) ! number of variables in node
      sz = scc(node)
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do piv = sptr(node), sptr(node+1)-1
         seen(piv) = node
         rlist(idx) = piv
         idx = idx + 1
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(j.lt.sptr(node)) cycle ! eliminated
            if(seen(j).eq.node) cycle ! already seen
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(j.lt.piv) cycle ! in upper triangle
            if(seen(j).eq.node) cycle ! already seen in this snode
            ! Otherwise, this is a new entry
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
      end do

      ! Note: following error check won't work with block pivots
      !if(idx .ne. rptr(node+1)) then
      !   ! Note: This bit of code should NEVER execute
      !  if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
      !      "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
      !      " entries, but expected to find ", rptr(node+1)-rptr(node)
      !   info = MC78_ERROR_UNKNOWN
      !   !print *, rlist(1:idx-1)
      !   return
      !endif
   end do
end subroutine mc78_row_lists_nosvar_integer

subroutine mc78_row_lists_svar_integer(nsvar, svar, n, ptr, row, perm, invp, &
      nnodes, sptr, sparent, scc, rptr, rlist, control, info, st)
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, intent(in) :: n
   integer(pkg_type), dimension(nsvar+1), intent(in) :: ptr
   integer, dimension(ptr(nsvar+1)-1), intent(in) :: row
   integer, dimension(nsvar), intent(in) :: perm
   integer, dimension(nsvar), intent(in) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer, dimension(nnodes), intent(in) :: scc
   integer(long), dimension(nnodes+1), intent(out) :: rptr
   integer, dimension(sum(scc(1:nnodes))), intent(out) :: rlist
   type(mc78_control), intent(in) :: control
   integer, intent(inout) :: info
   integer, intent(out) :: st

   integer :: child ! current child of node
   integer :: col ! current column of matrix corresponding to piv
   integer(long) :: i
   integer(long) :: idx ! current insert position into nodes(node)%index
   integer :: j
   integer :: k
   integer :: nelim ! number of variables eliminated at node
   integer :: node ! current node of assembly tree
   integer :: piv ! current pivot position
   integer, dimension(:), allocatable :: seen ! tag last time index was seen
   integer, dimension(:), allocatable :: chead ! head of child linked lists
   integer, dimension(:), allocatable :: cnext ! pointer to next child
   integer, dimension(:), allocatable :: svptr ! pointers for row list starts
   integer :: sz ! number of rows in node

   ! Allocate and initialise memory
   allocate(seen(n), chead(nnodes+1), cnext(nnodes+1), svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   seen(:) = 0
   chead(:) = -1

   ! Build svptr array
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(invp(i))
   end do

   ! Build child linked lists (backwards so pop off in good order)
   do node = nnodes, 1, -1
      i = sparent(node)
      cnext(node) = chead(i)
      chead(i) = node
   end do

   ! Loop over nodes from bottom up building row lists.
   rptr(1) = 1
   do node = 1, nnodes

      ! Allocate space for row indices
      nelim = sptr(node+1) - sptr(node) ! number of variables in node
      sz = scc(node)
      rptr(node+1) = rptr(node) + scc(node)
      idx = rptr(node) ! insert position

      ! Add entries eliminated at this node
      do i = sptr(node), sptr(node+1)-1
         do piv = svptr(i), svptr(i+1)-1
            seen(piv) = nnodes+1
            rlist(idx) = piv
            idx = idx + 1
         end do
      end do

      ! Find indices inherited from children
      child = chead(node)
      do while (child.ne.-1)
         do i = rptr(child), rptr(child+1)-1
            j = rlist(i)
            if(seen(j).ge.node) cycle ! already seen (or eliminated already)
            seen(j) = node
            rlist(idx) = j
            idx = idx + 1
         end do
         child = cnext(child)
      end do

      ! Find new indices from A
      do piv = sptr(node), sptr(node+1)-1
         col = invp(piv)
         do i = ptr(col), ptr(col+1)-1
            j = perm(row(i))
            if(seen(svptr(j)).ge.node) cycle ! already seen (or eliminated)
            ! Otherwise, this is a new entry
            ! Iterate over variables in supervariable
            do k = svptr(j), svptr(j+1)-1
               seen(k) = node
               rlist(idx) = k
               idx = idx + 1
            end do
         end do
      end do

      if(idx .ne. rptr(node+1)) then
         ! Note: This bit of code should NEVER execute
         if(control%unit_error.gt.0) write(control%unit_error, "(3(a,i8))") &
            "MC78 Internal Error: node ", node, ": found ", idx-rptr(node), &
            " entries, but expected to find ", rptr(node+1)-rptr(node)
         info = MC78_ERROR_UNKNOWN
         !print *, rlist(1:idx-1)
         return
      endif
   end do
end subroutine mc78_row_lists_svar_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Optimize cache locality routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! The following subroutine reorders the elimination order within each node in
! such a way that the order of the variables in
! the "primary" child has ordering agreement with the node. Consider
! the following tree:
!          A
!         / \
!        B   C
! 
! and assume B has more partially summed variables than C. Then
! B is the primary child and the ordering of the
! corresponding fully summed variables in the parent A matches the ordering
! of the partially summed variables in B (but not in C). Any additional
! partially summed variables present in C but not in B are then ordered in A
! such that they match C.
!
! This is done by two passes of the tree.
! The first builds a map from variables to the nodes at which
! they are eliminated, and orders the children of each node such
! that the first has the largest number of partially summed variables.
! The second pass uses a depth first search of the now ordered tree. It loops
! over non-fully summed variables and when it first encounters each it will
! place it as the next variable at its elimination node.
!
subroutine mc78_optimize_locality(n, realn, perm, invp, nnodes, sptr, sparent, &
      rptr, rlist, st, sort)

   integer, intent(in) :: n ! dimension of system
   integer, intent(in) :: realn ! symbolic dimension of system
   integer, dimension(n), intent(inout) :: perm ! on exit, will have been
      ! reordered for better cache locality
   integer, dimension(n), intent(inout) :: invp ! inverse of perm. on exit
      ! will have been changed to match new perm
   integer, intent(in) :: nnodes ! number of supernodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer, dimension(nnodes), intent(in) :: sparent
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st ! stat parameter
   logical, optional, intent(in) :: sort

   integer(long) :: i ! loop index
   integer :: id ! current insert position into list
   integer :: k
   integer :: j ! temporary variable
   integer, allocatable :: list(:) ! list of nodes in a weighted depth-first
      ! order such that children are visitied in decreasing number of
      ! partially summed variables (ie child with most p.s.v. visited first)
   integer, allocatable :: map(:) ! maps variables to nodes where
      ! they are eliminated
   integer :: node ! current node
   integer, allocatable :: ord(:) ! tracks number of variables ordered at
      ! each node
   integer, allocatable :: perm2(:) ! permutation to apply to perm
   integer :: pnode ! parent node
   integer :: shead ! current top of stack
   integer, allocatable :: stack(:) ! used for depth first walk of tree
   integer :: start ! first entry on stack of child from current node
   integer, allocatable :: chead(:) ! heads of child linked lists
   integer, allocatable :: cnext(:) ! tails of child linked lists

   ! Allocate arrays for depth first search of tree
   allocate (map(n), list(nnodes+1), stack(nnodes), chead(nnodes+1), &
      cnext(nnodes), stat=st)
   if (st /= 0) return

   !
   ! Build elimination map
   !
   do node = 1, nnodes
      do i = sptr(node), sptr(node+1)-1
         map(i) = node
      end do
   end do

   !
   ! Build child linked lists
   !
   chead(:) = -1 ! no child if necessary
   do i = nnodes, 1, -1 ! do in reverse order so they come off in original order
      j = sparent(i)
      cnext(i) = chead(j)
      chead(j) = i
   end do

   !
   ! Perform depth first search of tree, such that children of a node are
   ! visited in order of the number of partially summer variables, largest
   ! first.
   !
   shead = 1
   stack(shead) = nnodes + 1
   id = nnodes+1
   do while(shead.ne.0)
      ! Get node from top of stack
      node = stack(shead)
      shead = shead - 1

      ! Number it
      list(id) = node
      id = id - 1

      ! Place all its children on the stack
      start = shead + 1
      i = chead(node)
      do while(i.ne.-1)
         shead = shead + 1
         stack(shead) = i
         i = cnext(i)
      end do
      ! Order children just placed on stack such that child with least partially
      ! summed variables is at the top
      call order_children(shead-start+1, stack(start:shead), nnodes, sptr, &
         rptr, st)
      if(st.ne.0) return
   end do

   !
   ! Next loop over children reordering partially summed variables.
   !
   allocate(ord(nnodes),perm2(n),stat=st)
   if (st.ne.0) return

   do node = 1, nnodes
      ord(node) = sptr(node)
   end do

   do k = 1, nnodes
      node = list(k)

      ! Order variables first encountered at this node
      do i = rptr(node), rptr(node+1)-1
         j = rlist(i)
         pnode = map(j)
         if(pnode .ne. -1) then ! check if we have ordered j already
            ! order at parent
            perm2(j) = ord(pnode)
            ord(pnode) = ord(pnode) + 1
            map(j) = -1 ! mark as ordered
         endif
         rlist(i) = perm2(j)
      end do
   end do

   do i = realn+1, n
      perm2(i) = i
   end do

   !
   ! Apply permutation to perm and invp
   !
   ! Use perm as a temporary variable to permute invp.
   perm(1:n) = invp(1:n)
   do i = 1, n
      j = perm2(i)
      invp(j) = perm(i)
   end do

   ! Recover invp as inverse of perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   if(present(sort)) then
      if(sort) then
         call dbl_tr_sort(n, nnodes, rptr, rlist, st)
         if(st.ne.0) return
      endif
   endif
end subroutine mc78_optimize_locality

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Simple sort version, good for nodes with small numbers of children
! (Passes to mergesort for large numbers of entries)
recursive subroutine order_children(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: nelim, m
   integer :: k, kdummy, klo, kor
   integer :: ice_idx, ice_psum, ik_idx, ik_psum

   st = 0

   if(n.ge.minsz_ms) then
      call order_children_ms(n, child, nnodes, sptr, rptr, st)
   else
      klo = 2
      kor = n
      do kdummy = klo, n
         ! items kor, kor+1, .... ,n are in order
         ice_idx = child(kor-1)
         nelim = sptr(ice_idx+1) - sptr(ice_idx)
         m = rptr(ice_idx+1) - rptr(ice_idx)
         ice_psum = m - nelim
         do k = kor, n
            ik_idx = child(k)
            nelim = sptr(ik_idx+1) - sptr(ik_idx)
            m = rptr(ik_idx+1) - rptr(ik_idx)
            ik_psum = m - nelim
            if (ice_psum .ge. ik_psum) exit
            child(k-1) = ik_idx
         end do
         child(k-1) = ice_idx
         kor = kor - 1
      end do
   endif
end subroutine order_children

! Orders nodes stored in child(1:n) such that the number of partially summed
! variables at each node is decreasing (ie one with least is in posn n)
!
! Merge sort version, dramatically improves performance for nodes with large
! numbers of children
! (Passes to simple sort for small numbers of entries)
recursive subroutine order_children_ms(n, child, nnodes, sptr, rptr, st)
   integer, intent(in) :: n
   integer, dimension(n), intent(inout) :: child
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, intent(out) :: st

   integer :: i, j, jj, jj2, k, kk, kk2, m, nelim
   integer :: mid
   integer, dimension(:), allocatable :: work

   if(n.le.1) return
   if(n.lt.minsz_ms) then
      call order_children(n, child, nnodes, sptr, rptr, st)
      return
   endif
   mid = (n-1)/2 + 1

   ! Recurse to order half lists
   call order_children_ms(mid, child(1:mid), nnodes, sptr, rptr, st)
   if(st.ne.0) return
   call order_children_ms(n - mid, child(mid+1:n), nnodes, sptr, rptr, st)
   if(st.ne.0) return

   ! Merge two half lists
   ! (Take a copy of the first half list so we don't overwrite it)
   allocate(work(mid), stat=st)
   if(st.ne.0) return
   work(:) = child(1:mid)
   j = 1
   k = mid+1
   jj = work(j)
   nelim = sptr(jj+1) - sptr(jj)
   m = rptr(jj+1) - rptr(jj)
   jj2 = m - nelim
   kk = child(k)
   nelim = sptr(kk+1) - sptr(kk)
   m = rptr(kk+1) - rptr(kk)
   kk2 = m - nelim
   do i = 1, n
      if(jj2.ge.kk2) then
         child(i) = jj
         j = j + 1
         if(j.gt.mid) exit
         jj = work(j)
         nelim = sptr(jj+1) - sptr(jj)
         m = rptr(jj+1) - rptr(jj)
         jj2 = m - nelim
      else
         child(i) = kk
         k = k + 1
         if(k.gt.n) exit
         kk = child(k)
         nelim = sptr(kk+1) - sptr(kk)
         m = rptr(kk+1) - rptr(kk)
         kk2 = m - nelim
      endif
   end do
   if(j.le.mid) child(i+1:n) = work(j:mid)
end subroutine order_children_ms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assorted auxilary routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Converts piv_size to block_pivots:
! piv_size(i) is size of block pivot containing column i of A
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
!
subroutine convert_to_blk_piv(n, invp, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: invp
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, cnt

   allocate(blk2(n))

   ! Take a copy of block so we can change it
   blk2(1:n) = block(1:n)

   ! Iterate over pivots in elimination order, recording starts and ends of blks
   cnt = blk2(invp(1))-1 ! Initialise for first pivot
   block(1) = 1 ! First pivot is start of a block
   do i = 2, n
      block(i) = 0
      if(cnt.eq.0) then
         ! this is first pivot of a block, previous is last pivot of a block
         cnt = blk2(invp(i))
         block(i-1) = block(i-1) + 2
         block(i) = block(i) + 1
      endif
      cnt = cnt - 1
   end do
   block(n) = block(n) + 2 ! end of matrix must end a block pivot

end subroutine convert_to_blk_piv

!
! Converts block_pivots back to piv_size:
! block_pivots(j) is a flag for pivot j of L and has one of the following values
!  0 - pivot i is in the middle of a block pivot
!  1 - pivot i is the first pivot of a block pivot
!  2 - pivot i is the last pivot of a block pivot
!  3 - pivot i is a 1x1 pivot
! piv_size(i) is size of block pivot containing column i of A
!
subroutine convert_from_blk_piv(n, perm, block)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: block

   integer, dimension(:), allocatable :: blk2
   integer :: i, sa, cnt

   allocate(blk2(n))

   ! convert first/last notation to block size notation
   cnt = -1; sa = -1 ! these values should never actually be used
   do i = 1, n
      select case(block(i))
      case (0) ! middle pivot of a block
         cnt  = cnt + 1
      case (1) ! first pivot of a block
         sa = i
         cnt = 1
      case (2) ! end pivot of a block
         cnt = cnt + 1
         block(sa:i) = cnt
      case (3) ! only pivot of a block
         block(i) = 1
      end select
   end do

   ! Permute back to original matrix order
   blk2(1:n) = block(1:n)
   do i = 1, n
      block(i) = blk2(perm(i))
   end do
end subroutine convert_from_blk_piv

!
! This subroutine copies a matrix pattern and adds subdiagonal entries as
! needed to force block pivots to have a parent-child relation in the
! elimination tree.
!
subroutine mc78_block_prep(n, ptr, row, bptr, brow, perm, invp, block_pivots, &
      st)
   integer, intent(in) :: n ! Dimension of system
   integer(pkg_type), dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer(pkg_type), dimension(n+1), intent(out) :: bptr ! Column pointers
   integer, dimension(:), intent(out) :: brow ! Row indices
   integer, dimension(n), intent(inout) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
      ! On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(inout) :: invp ! inverse permutation of perm
   integer, dimension(n), intent(inout) :: block_pivots ! Matches pivot order
      ! and specifies block pivots.
      ! block_pivots(i) may take one of the following values:
      ! 0 - pivot i is in the middle of a block pivot
      ! 1 - pivot i is the first pivot of a block pivot
      ! 2 - pivot i is the last pivot of a block pivot
      ! 3 - pivot i is a 1x1 pivot
   integer, intent(out) :: st

   integer :: col
   integer :: i
   integer :: idx
   integer(pkg_type) :: j
   integer :: k
   integer :: piv
   integer, dimension(:), allocatable :: seen

   allocate(seen(n), stat=st)
   if(st.ne.0) return

   ! First pass through ptr and ensure that all block pivots contain no
   ! empty columns
   perm(:) = invp(:)
   piv = 1
   j = 1
   ! Handle variables that are actually used
   do while(piv.le.n)
      do i = piv, n
         if(block_pivots(i).ge.2) exit ! end of block pivot
      end do

      k = 0
      do piv = piv, i
         if(ptr(perm(piv)).eq.ptr(perm(piv)+1)) cycle
         invp(j) = perm(piv)
         j = j + 1
         if(k.eq.0) then
            ! This is the new start of the block pivot
            select case(block_pivots(piv))
            case(0) ! was in the middle. now a start
               block_pivots(piv) = 1
            case(2) ! was the end. now a 1x1
               block_pivots(piv) = 3
            end select
         endif
         k = piv
      end do
      if(k.ne.0) then
         ! The was at least one used variable in the block pivot
         select case(block_pivots(k))
         case(0) ! was the middle. now an end
            block_pivots(k) = 2
         case(1) ! was the start. now a 1x1
            block_pivots(k) = 3
         end select
      endif
      piv = i + 1
   end do
   ! Handle unused variables
   do piv = 1, n
      i = perm(piv)
      if(ptr(i).eq.ptr(i+1)) then
         invp(j) = i
         j = j + 1
         block_pivots(piv) = 3 ! Force to 1x1
      endif
   end do
   ! Map block_pivots in original variable order into sv_map
   do i = 1, n
      seen(perm(i)) = block_pivots(i)
   end do
   ! Map sv_map in new pivot order back into block_pivots
   do i = 1, n
      block_pivots(i) = seen(invp(i))
   end do
   ! Reestablish perm
   do i = 1, n
      perm(invp(i)) = i
   end do

   ! Now iterate over cleaned up block pivot sequence
   seen(:) = 0
   piv = 1
   idx = 1
   do col = 1, n
      piv = perm(col)
      bptr(col) = idx
      if(block_pivots(piv).eq.3) then
         ! 1x1 pivot, just copy the column
         idx = idx + ptr(col+1) - ptr(col)
         brow(bptr(col):idx-1) = row(ptr(col):ptr(col+1)-1)
      else
         ! copy the column, but add an entry on subdiagonal(s)
         do i = ptr(col), ptr(col+1)-1
            j = row(i)
            seen(j) = col
            brow(idx) = j
            idx = idx + 1
         end do
         if(block_pivots(piv).ne.1) then
            ! Not the first column, add an entry above the diagonal
            j = invp(piv-1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
         if(block_pivots(piv).ne.2) then
            ! Not the last column, add an entry below the diagonal
            j = invp(piv+1)
            if(seen(j).lt.col) then
               brow(idx) = j
               idx = idx + 1
            endif
         endif
      endif
   end do
   bptr(n+1) = idx
end subroutine mc78_block_prep

!
! This subroutine will take information concerning a compressed matrix and a
! supervariable map, and will decompress the information so it relates to
! the original matrix
!
subroutine svar_unmap(n, nsvar, svar, perm, invp, nnodes, sinvp, &
      snptr, st)
   integer, intent(in) :: n
   integer, intent(in) :: nsvar
   integer, dimension(nsvar), intent(in) :: svar
   integer, dimension(n), intent(out) :: perm
   integer, dimension(n), intent(inout) :: invp
   integer, intent(in) :: nnodes
   integer, dimension(nsvar), intent(in) :: sinvp
   integer, dimension(nnodes+1), intent(inout) :: snptr
   integer, intent(out) :: st

   integer, dimension(:), allocatable :: svptr
   integer :: i, j, k
   integer :: j1, j2
   integer :: idx

   ! Set up svptr
   allocate(svptr(nsvar+1), stat=st)
   if(st.ne.0) return
   svptr(1) = 1
   do i = 1, nsvar
      svptr(i+1) = svptr(i) + svar(i)
   end do

   ! Take a copy of invp in perm to ease remapping
   perm(:) = invp(:)

   ! Remap invp
   idx = 1
   do i = 1, nsvar
      j = sinvp(i)
      do k = svptr(j), svptr(j+1)-1
         invp(idx) = perm(k)
         idx = idx + 1
      end do
   end do

   ! Expand supernode pointer
   j1 = snptr(1)
   do i = 1, nnodes
      j2 = snptr(i+1)
      snptr(i+1) = snptr(i)
      do j = j1, j2-1
         snptr(i+1) = snptr(i+1) + svar(sinvp(j))
      end do
      j1 = j2
   end do

   ! Finally, recover perm as inverse of invp
   do i = 1, n
      perm(invp(i)) = i
   end do
end subroutine svar_unmap


!
! This subroutine performs a double transpose sort on the row indices of sn
!
subroutine dbl_tr_sort(n, nnodes, rptr, rlist, st)
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(inout) :: rlist
   integer, intent(out) :: st

   integer :: node
   integer(long) :: i
   integer(long) :: j
   integer(long), dimension(:), allocatable :: ptr
   integer, dimension(:), allocatable :: nptr
   integer, dimension(:), allocatable :: col

   allocate(ptr(n+2), stat=st)
   if(st.ne.0) return
   ptr(:) = 0

   ! Count number of entries in each row. ptr(i+2) = #entries in row i
   do node = 1, nnodes
      do i = rptr(node), rptr(node+1)-1
         j = rlist(i) ! row entry
         ptr(j+2) = ptr(j+2) + 1
      end do
   end do

   ! Determine row starts. ptr(i+1) = start of row i
   ptr(1:2) = 1
   do i = 1, n
      ptr(i+2) = ptr(i+1) + ptr(i+2)
   end do

   j = ptr(n+2)-1 ! total number of entries
   allocate(col(j), stat=st)
   if(st.ne.0) return

   ! Now fill in col array
   do node = 1, nnodes
      do i = rptr(node), rptr(node+1)-1
         j = rlist(i) ! row entry
         col( ptr(j+1) ) = node
         ptr(j+1) = ptr(j+1) + 1
      end do
   end do

   ! Finally transpose back into nodes
   allocate(nptr(nnodes))
   nptr(:) = rptr(1:nnodes)
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         node = col(j)
         rlist(nptr(node)) = i
         nptr(node) = nptr(node) + 1
      end do
   end do
end subroutine dbl_tr_sort

!
! This subroutine applies the permutation perm to order, invp and cc
!
subroutine apply_perm(n, perm, order, invp, cc, block_pivots)
   integer, intent(in) :: n
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(inout) :: order
   integer, dimension(n), intent(inout) :: invp
   integer, dimension(n), intent(inout) :: cc
   integer, dimension(n), optional, intent(inout) :: block_pivots

   integer :: i
   integer :: j

   ! Use order as a temporary variable to permute cc. Don't care about cc(n+1)
   order(1:n) = cc(1:n)
   do i = 1, n
      j = perm(i)
      cc(j) = order(i)
   end do

   ! Use order as a temporary variable to permute invp.
   order(1:n) = invp(1:n)
   do i = 1, n
      j = perm(i)
      invp(j) = order(i)
   end do

   ! Use order as a temporary variable to permute block_pivots if present
   if(present(block_pivots)) then
      order(1:n) = block_pivots(1:n)
      do i = 1, n
         j = perm(i)
         block_pivots(j) = order(i)
      end do
   endif

   ! Recover order as inverse of invp
   do i = 1, n
      order(invp(i)) = i
   end do
end subroutine apply_perm

end module hsl_mc78_integer
! COPYRIGHT (c) 2009 Council for the Central Laboratory
!               of the Research Councils
! Original date 20 October 2009. Version 1.0.0.

! Fortran 95 version of the mc34 package.

! 18 May 2010 Version 1.1.0 -jhogg
!             Create hsl_mc34_integer
!             Change from logical Hermitian to integer sym_type to cope with
!             skew symmetric matrices as well as Hermitian and symmetric ones.

! to change precision:
!    change _double, kind(0.0d0)
! For complex version:
!    change real to complex
! For Hermitian case, take conjugate for upper triangle entries
! For integer version:
!    change real to integer, kind(0)

   module hsl_mc34_double
   implicit none
   private
   public mc34_expand

   integer, parameter :: wp = kind(0.0d0)

   interface mc34_expand
      module procedure mc34_expand_double
   end interface

   contains

      subroutine mc34_expand_double(n,row,ptr,iw,a,sym_type)

!  this subroutine generates the expanded structure for a
!   matrix a with a symmetric sparsity pattern given the structure 
!   for the lower triangular part.  diagonal entries need not be present.

      integer, intent(in) :: n  ! holds the order of a.

      integer, intent(inout) :: row(*) ! must be set by the user to
!       hold the row indices of the lower triangular part of a.
!       the entries of a single column must be
!       contiguous. the entries of column j must precede those of column
!       j+1, and there must be no wasted space between
!       columns. row indices within a column may be in any order.  on
!       exit, it will have the same meaning but will be changed to hold
!       the row indices of the entries in the expanded structure.  diagonal
!       entries need not be present. the new row indices added in the
!       upper triangular part will be in order for each column and will
!       precede the row indices for the lower triangular part which will
!       remain in the input order.

     integer, intent(inout) ::ptr(n+1)  !  must be set
!       by the user so that ptr(j) is the position in row
!       of the first entry in column j and
!       ptr(n+1) must be set to one more than the total number of
!       entries.  on exit, ptr(j) will have the same meaning but
!       will be changed to point to the position of the first entry of
!       column j in the expanded structure. the new value of
!       ptr(n+1) will be one greater than the number of entries in
!       the expanded structure.

     integer :: iw(n) ! workspace

     real(wp), optional, intent(inout) :: a(*) 
!       if present, a(1:ptr(n+1)-1) must be set by the user so that
!       a(k) holds the value of the entry in row(k). 
!       on exit, a will hold the values of the entries in the expanded 
!       structure corresponding to the output values of row.

     integer, optional, intent(in) :: sym_type 
!      if present with value 1, matrix is skew symmetric.
!      if present with value 2, matrix is hermitian.
!      otherwise matrix is symmetric.

      integer :: ckp1 ! used as running pointer
      integer :: i,i1,i2,ii,ipkp1,ipos
      integer :: j,jstart 
      integer :: lenk ! number of entries in col. j of original structure
      integer :: ndiag ! number diagonal entries present
      integer :: newtau ! number of entries in expanded storage
      integer :: oldtau ! number of entries in symmetric storage
      integer :: r_sym_type ! real sym_type value (used as argument is optional)

      oldtau = ptr(n+1) - 1
      iw(1:n) = 0

! iw(j) set to total number entries in col. j of expanded mx.
      ndiag = 0
      do j = 1,n
        i1 = ptr(j)
        i2 = ptr(j+1) - 1
        iw(j) = iw(j) + i2 - i1 + 1
        do ii = i1,i2
          i = row(ii)
          if (i /= j) then
            iw(i) = iw(i) + 1
          else
            ndiag = ndiag + 1
          end if
        end do
      end do

      newtau = 2*oldtau - ndiag
! ipkp1 points to position  after end of column being currently processed
      ipkp1 = oldtau + 1
! ckp1 points to position  after end of same column in expanded structure
      ckp1 = newtau + 1
! go through the array in the reverse order placing lower triangular
!     elements in  appropriate slots.
      do j = n,1,-1
        i1 = ptr(j)
        i2 = ipkp1
        lenk = i2 - i1
! jstart is running pointer to position in new structure
        jstart = ckp1
! set ikp1 for next column
        ipkp1 = i1
        i2 = i2 - 1
! run through columns in reverse order
! lower triangular part of col. moved to end of same column in expanded form
        if (present(a)) then
          do ii = i2,i1,-1
            jstart = jstart - 1
            a(jstart) = a(ii)
            row(jstart) = row(ii)
          end do
        else
          do ii = i2,i1,-1
            jstart = jstart - 1
            row(jstart) = row(ii)
          end do
        end if
! ptr is set to position of first entry in lower triangular part of
!     column j in expanded form
        ptr(j) = jstart
! set ckp1 for next column
        ckp1 = ckp1 - iw(j)
! reset iw(j) to number of entries in lower triangle of column.
        iw(j) = lenk
      end do
 
! again sweep through the columns in the reverse order, this
!     time when one is handling column j the upper triangular
!     elements a(j,i) are put in position.
        do j = n,1,-1
          i1 = ptr(j)
          i2 = ptr(j) + iw(j) - 1
! run down column in order
! note that i is always greater than or equal to j
          if (present(a)) then
            r_sym_type = 0 ! symmetric
            if(present(sym_type)) r_sym_type = sym_type
            select case(r_sym_type)
            case(1) ! skew symmetric
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = -a(ii)
                row(ipos) = j
              end do
            case default ! symmetric or hermitian
              do ii = i1,i2
                i = row(ii)
                if (i == j) cycle
                ptr(i) = ptr(i) - 1
                ipos = ptr(i)
                a(ipos) = a(ii)
                row(ipos) = j
              end do
            end select
          else
            do ii = i1,i2
              i = row(ii)
              if (i == j) cycle
              ptr(i) = ptr(i) - 1
              ipos = ptr(i)
              row(ipos) = j
            end do
          end if
        end do
      ptr(n+1) = newtau + 1

      end subroutine mc34_expand_double
   end module hsl_mc34_double
! COPYRIGHT (c) 2007 Science & Technology Facilities Council
! Original date 1 August 2007. Version 1.0.0. (Sue Dollar)
!
! 1 Decemeber 2010 Version 2.0.0 (Jonathan Hogg)
!    Modify interface to allow specificaiton of source and destination ranges,
!    substantially rewrite and simplify code, add long integer support, remove
!    support for unallocated arrays on input

! To convert to single:
!    s/double/single
!    change myreal definition
! To convert to integer:
!    s/double/integer
!    s/real(myreal)/integer(myinteger)
!    s/REAL (kind=myreal)/INTEGER (kind=myinteger)
!    change myreal definintion to myinteger
! To convert to long integer:
!    s/double/long
!    s/real(myreal)/integer(myinteger)
!    s/REAL (kind=myreal)/INTEGER (kind=myinteger)
!    change myreal definintion to myinteger
! To convert to double complex:
!    s/double/complex
!    s/COMPLEX (kind=mycomplex)/COMPLEX (kind=mycomplex)
!    s/complex(mycomplex)/complex(mycomplex)
!    change myreal definition to mycomplex
! To convert to complex:
!    s/double/complex
!    s/COMPLEX (kind=mycomplex)/COMPLEX (kind=mycomplex)
!    s/complex(mycomplex)/complex(mycomplex)
!    change myreal definition to mycomplex
MODULE hsl_zb01_integer

   IMPLICIT NONE
   PRIVATE

   ! ---------------------------------------------------
   ! Precision
   ! ---------------------------------------------------

   INTEGER, PARAMETER :: myinteger = kind(1)
   INTEGER, PARAMETER :: myint = kind(1)
   INTEGER, PARAMETER :: long = selected_int_kind(18)

   ! ---------------------------------------------------
   ! Error flags
   ! ---------------------------------------------------
   INTEGER (kind=myint), PARAMETER :: &
      zb01_err_lw = -1, &        ! lw<=lkeep/size(w) on input
      zb01_err_lw_1 = -2, &      ! lw<1 on input
      zb01_err_lw_both = -3, &   ! both -1 and -2
      zb01_err_lkeep = -4, &     ! lkeep>size(w) on input
      zb01_err_lkeep_l = -5, &   ! lkeep<1
      zb01_err_lkeep_both = -6, &! both -4 and -5
      zb01_err_filename = -7, &  ! filename too long
      zb01_err_file_size = -8, & ! file_size <2**12
      zb01_err_filename_exists = -9, & ! filename already exists
      zb01_err_mode = -10, &     ! mode out of range
      zb01_err_memory_alloc = -11, & ! memory alloc error
      zb01_err_memory_dealloc = -12, & ! memory dealloc error
      zb01_err_inquire = -13, &  ! error in Fortran inquire statement
      zb01_err_open = -14, &     ! error in Fortran open statement
      zb01_err_read = -15, &     ! error in Fortran read statement
      zb01_err_write = -16, &    ! error in Fortran write statement
      zb01_err_close = -17, &    ! error in Fortran close statement
      zb01_err_src_dest = -18, & ! src and dest sizes do not match
      zb01_err_w_unalloc = -19   ! w is unallocated on entry

   ! ---------------------------------------------------
   ! Warning flags
   ! ---------------------------------------------------
   INTEGER (kind=myint), PARAMETER :: &
      zb01_warn_lw = 1 ! size_out changed

   ! ---------------------------------------------------
   ! Derived type definitions
   ! ---------------------------------------------------

   TYPE, PUBLIC :: zb01_info
      INTEGER :: flag = 0        ! error/warning flag
      INTEGER :: iostat = 0      ! holds Fortran iostat parameter
      INTEGER :: stat = 0        ! holds Fortran stat parameter
      INTEGER :: files_used = 0  ! unit number scratch file written to
   END TYPE zb01_info

   INTERFACE zb01_resize1
      MODULE PROCEDURE zb01_resize1_integer
   END INTERFACE

   INTERFACE zb01_resize2
      MODULE PROCEDURE zb01_resize2_integer
   END INTERFACE

   PUBLIC zb01_resize1, zb01_resize2

CONTAINS

SUBROUTINE zb01_resize1_integer(w,size_in,size_out,info,src,dest,filename, &
      file_size,mode)

  ! --------------------------------------
  ! Expands w to have larger size and copies all or part of the
  ! original array into the new array
  ! --------------------------------------

  ! w: is a REAL allocatable array of INTENT(INOUT). Must be allocated on entry.
  INTEGER (kind=myinteger), DIMENSION (:), ALLOCATABLE, INTENT (INOUT) :: w

  ! size_in: is an INTEGER OF INTENT(IN). It holds the extent of w on entry.
  INTEGER (kind=long), INTENT(IN) :: size_in

  ! size_out: is an INTEGER of INTENT(INOUT). It holds the required size of w
  ! on input and the actual size of w on output
  INTEGER (kind=long), INTENT (INOUT) :: size_out

  ! info: is of derived type ZB01_info with intent(out).
  ! info%flag = ZB01_ERR_LW if 1<=size_out<=lkeep/size(w) on input
  ! ZB01_ERR_LW_L if size_out<=0
  ! ZB01_ERR_LKEEP if lkeep>size(w) on input
  ! ZB01_ERR_LKEEP_L if lkeep<1 on input
  ! ZB01_ERR_FILENAME if filename too long
  ! ZB01_ERR_FILE_SIZE if file_size <2**12
  ! ZB01_ERR_FILENAME_EXISTS if filename already exists
  ! ZB01_ERR_MODE if mode out of range
  ! ZB01_ERR_MEMORY_ALLOC if memory alloc error
  ! ZB01_ERR_MEMORY_DEALLOC if memory dealloc error
  ! ZB01_ERR_INQUIRE if error in Fortran inquire statement
  ! ZB01_ERR_OPEN if error in Fortran open statement
  ! ZB01_ERR_READ if error in Fortran read statement
  ! ZB01_ERR_WRITE if error in Fortran write statement
  ! ZB01_ERR_CLOSE if error in Fortran close statement
  ! ZB01_WARN_LW if size_out changed by subroutine
  ! ZB01_WARN_W if w not allocated on input
  ! info%iostat holds Fortran iostat parameter
  ! info%stat holds Fortran stat parameter
  ! info%unit holds unit number scratch file written to (negative if
  ! not)
  TYPE (zb01_info), INTENT (OUT) :: info

  ! src and dest: are OPTIONAL INTEGER arrays of INTENT(IN). They specify
  ! the source and destination ranges for any values to be kept.
  ! dest(2)-dest(1) must equal src(2)-src(1) and minval(src,dest)>0.
  INTEGER (kind=long), DIMENSION(2), INTENT (IN), OPTIONAL :: src
  INTEGER (kind=long), DIMENSION(2), INTENT (IN), OPTIONAL :: dest

  ! filename: is an OPTIONAL STRING OF CHARACTERS of INTENT(IN).
  ! It holds the name of the file to be used
  CHARACTER (len=*), INTENT (IN), OPTIONAL :: filename

  ! file_size: is an OPTIONAL INTEGER of INTENT(IN). It holds the length
  ! of
  ! of the files to be used if necessary
  INTEGER (kind=long), INTENT (IN), OPTIONAL :: file_size

  ! mode: is an OPTIONAL INTEGER of INTENT(IN). If mode==0, then the
  ! subroutine will firstly try to use a temporary array. If mode==1,
  ! then
  ! the subroutine will immediately use temporary files
  INTEGER, INTENT (IN), OPTIONAL :: mode

  ! --------------------------------------
  ! Local variables
  ! --------------------------------------

  ! wtemp: is a REAL allocatable array
  INTEGER (kind=myinteger), DIMENSION (:), ALLOCATABLE :: wtemp

  ! length to use for wtemp
  INTEGER :: lwtemp

  ! Fortran stat parameter
  INTEGER :: stat

  INTEGER :: number_files, mode_copy
  INTEGER (kind=long) :: file_size_copy
  INTEGER, DIMENSION (:), ALLOCATABLE :: units
  integer (kind=long), dimension(2,2) :: src_copy, dest_copy

  ! --------------------------------------
  ! Check input for errors
  ! --------------------------------------

  info%flag = 0

  ! Check w is allocated
  if(.not.allocated(w)) then
    info%flag = zb01_err_w_unalloc
    return
  endif

  ! Setup src_copy and dest_copy
  src_copy(1,1) = 1
  src_copy(2,1) = size_in
  src_copy(1,2) = 1
  src_copy(2,2) = 1
  if(present(src)) then
    src_copy(:,1) = src(:)
  elseif(present(dest)) then
    src_copy(:,1) = dest(:)
  endif
  dest_copy(:,:) = src_copy(:,:)
  if(present(dest)) then
    dest_copy(:,1) = dest(:)
  endif

  if(minval(src_copy).lt.0 .or. minval(dest_copy).lt.0) &
    info%flag = zb01_err_lkeep_l

  if(maxval(src_copy) .gt. size_in) then
    if(info%flag.eq.zb01_err_lkeep_l) then
      info%flag = zb01_err_lkeep_both
    else
      info%flag = zb01_err_lkeep
    endif
  endif
  if(info%flag.lt.0) return
  
  if(any(src_copy(2,:)-src_copy(1,:) .ne. &
      dest_copy(2,:)-dest_copy(1,:))) then
    info%flag = zb01_err_src_dest
    return
  endif

  if(size_out.lt.1) then
    info%flag = zb01_err_lw_1
    return
  endif
  if(maxval(dest_copy).gt.size_out) then
    info%flag = zb01_err_lw
    return
  endif

  IF (present(filename)) THEN
    IF (len(filename)>400) THEN
      ! Error: len(filename) > 400
      info%flag = zb01_err_filename
      RETURN
    END IF
  END IF


  file_size_copy = 2**22_long
  if (present(file_size)) file_size_copy = file_size

  IF (file_size_copy<2**12) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_file_size
    RETURN
  END IF

  mode_copy = 0
  if(present(mode)) mode_copy = mode
  IF (mode_copy<0 .OR. mode_copy>1) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_mode
    RETURN
  END IF

  IF (src_copy(2,1)-src_copy(1,1)+1.le.0) THEN

    ! --------------------------------------
    ! No entries need to be copied
    ! --------------------------------------

    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out==size_in) RETURN

    ! --------------------------------------
    ! Deallocate w
    ! --------------------------------------
    DEALLOCATE (w,STAT=info%stat)
    IF (info%stat>0) THEN
      info%flag = zb01_err_memory_dealloc
      RETURN
    END IF

    ! --------------------------------------
    ! Reallocate w
    ! --------------------------------------
    ALLOCATE (w(size_out),STAT=info%stat)
    IF (info%stat>0) THEN
      ! --------------------------------------
      ! Allocation of w failed
      ! --------------------------------------
      info%flag = zb01_err_memory_alloc
      RETURN
    END IF

  ELSE

    ! --------------------------------------
    ! entries need to be copied
    ! --------------------------------------


    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out==size_in .and. all(src_copy(:,:).eq.dest_copy(:,:))) RETURN

    ! --------------------------------------
    ! Attempt to allocate temporary array
    ! --------------------------------------
    stat = 0
    ! Length of temporary array
    lwtemp = src_copy(2,1) - src_copy(1,1) + 1

    ! Allocate temporary array
    IF (mode_copy==0) ALLOCATE (wtemp(lwtemp),STAT=stat)

    IF (stat==0 .AND. mode_copy==0) THEN
      ! --------------------------------------
      ! Allocation successful
      ! --------------------------------------

      ! --------------------------------------
      ! Copy entries into wtemp
      ! --------------------------------------
      wtemp(1:lwtemp) = w(src_copy(1,1):src_copy(2,1))
      !print *, "copy w to wtemp"

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out),STAT=stat)
      IF (stat>0) THEN
        !print *, "fail alloc of w, copy wtemp to file", src_copy, size_out


        ! --------------------------------------
        ! Allocation not successful
        ! --------------------------------------
        call write_to_file(wtemp, size_in, src_copy, &
          units, number_files, file_size_copy, info, filename)
        if(info%flag.lt.0) return

        ! --------------------------------------
        ! Deallocate wtemp
        ! --------------------------------------
        DEALLOCATE (wtemp,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = zb01_err_memory_dealloc
          RETURN
        END IF

        ! --------------------------------------
        ! Reallocate w
        ! --------------------------------------
        ALLOCATE (w(size_out),STAT=info%stat)
        IF (info%stat>0) THEN
          !print *, "still can't alloc w. give up (try and preserve)."
          ! --------------------------------------
          ! Allocation of w to desired size failed
          ! Try size=lwtemp
          ! --------------------------------------


          ALLOCATE (w(lwtemp),STAT=info%stat)
          IF (info%stat==0) THEN
            ! --------------------------------------
            ! Reassign lw
            ! --------------------------------------
            size_out = lwtemp
            info%flag = zb01_warn_lw

          ELSE
            !print *, "can't even do that!"

            ! --------------------------------------
            ! Allocation of w still failed - delete files
            ! --------------------------------------
            call delete_files(units, number_files, info, filename)
            if(info%flag.lt.0) return

            info%flag = zb01_err_memory_alloc
            RETURN
          END IF
        END IF

        !print *, "restore entries to w"
        ! --------------------------------------
        ! Copy entries back into w
        ! --------------------------------------
        call read_from_file(w, size_in, dest_copy, units, &
          number_files, file_size_copy, info, filename)

        RETURN
      END IF

      !print *, "copy wtemp back to w"
      ! --------------------------------------
      ! Copy entries
      ! --------------------------------------
      w(dest_copy(1,1):dest_copy(2,1)) = wtemp(1:lwtemp)

      ! --------------------------------------
      ! Deallocate wtemp
      ! --------------------------------------
      DEALLOCATE (wtemp,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF
      RETURN

    ELSE
      !print *, "can't alloc wtemp, stow direct in file"
      ! --------------------------------------
      ! Allocation not successful
      ! --------------------------------------
      call write_to_file(w, size_in, src_copy, units, &
        number_files, file_size_copy, info, filename)
      if(info%flag.lt.0) return

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out),STAT=info%stat)
      IF (info%stat>0) THEN
        !print *, "can't alloc w"
        ! --------------------------------------
        ! Allocation of w to desired size failed
        ! Try size=lwtemp
        ! --------------------------------------

        ALLOCATE (w(lwtemp),STAT=info%stat)
        IF (info%stat==0) THEN
          ! --------------------------------------
          ! Reassign lw
          ! --------------------------------------
          size_out = lwtemp
          info%flag = zb01_warn_lw

        ELSE
          ! --------------------------------------
          ! Allocation of w failed again - delete files
          ! --------------------------------------
          call delete_files(units, number_files, info, filename)
          if(info%flag.lt.0) return

          info%flag = zb01_err_memory_alloc
          RETURN
        END IF
      END IF

      !print *, "restore values from file"
      call read_from_file(w, size_in, dest_copy, units, &
        number_files, file_size_copy, info, filename)

    END IF
  END IF

END SUBROUTINE zb01_resize1_integer

! ---------------------------------------------------------------

SUBROUTINE zb01_resize2_integer(w,size_in,size_out,info,src,dest,filename, &
      file_size,mode)

  ! --------------------------------------
  ! Expands w to have larger size and copies all or part of the
  ! original array into the new array
  ! --------------------------------------

  ! w: is a REAL allocatable array of INTENT(INOUT). Must be allocated on entry.
  INTEGER (kind=myinteger), DIMENSION (:,:), ALLOCATABLE, INTENT (INOUT) :: w

  ! size_in: is an INTEGER array INTENT(IN). It holds the extent of w on entry.
  INTEGER (kind=long), INTENT(IN) :: size_in(2)

  ! size_out: is an INTEGER array of rank-one with size 2 and of
  ! INTENT(INOUT).
  ! On input, size_out(1) holds the required number of rows of w and size_out(2)
  ! holds
  ! the required number of columns of w. On successful output it
  ! contains
  ! the number of rows and columns that w has.
  INTEGER (kind=long), INTENT (INOUT) :: size_out(2)

  ! info: is of derived type ZB01_info with intent(out).
  ! info%flag = ZB01_ERR_LW if 1<=size_out<=lkeep/size(w) on input
  ! ZB01_ERR_LW_L if size_out<=0
  ! ZB01_ERR_LKEEP if lkeep>size(w) on input
  ! ZB01_ERR_LKEEP_L if lkeep<1 on input
  ! ZB01_ERR_FILENAME if filename too long
  ! ZB01_ERR_FILE_SIZE if file_size <2**12
  ! ZB01_ERR_FILENAME_EXISTS if filename already exists
  ! ZB01_ERR_MODE if mode out of range
  ! ZB01_ERR_MEMORY_ALLOC if memory alloc error
  ! ZB01_ERR_MEMORY_DEALLOC if memory dealloc error
  ! ZB01_ERR_INQUIRE if error in Fortran inquire statement
  ! ZB01_ERR_OPEN if error in Fortran open statement
  ! ZB01_ERR_READ if error in Fortran read statement
  ! ZB01_ERR_WRITE if error in Fortran write statement
  ! ZB01_ERR_CLOSE if error in Fortran close statement
  ! ZB01_WARN_LW if size_out changed by subroutine
  ! ZB01_WARN_W if w not allocated on input
  ! info%iostat holds Fortran iostat parameter
  ! info%stat holds Fortran stat parameter
  ! info%files_used holds unit number of files written to
  TYPE (zb01_info), INTENT (OUT) :: info

  ! src and dest: are OPTIONAL INTEGER arrays of INTENT(IN). They specify
  ! the source and destination ranges for any values to be kept.
  ! dest(2,:)-dest(1,:) must equal src(2,:)-src(1,:) and must be
  ! non-negative.
  INTEGER (kind=long), DIMENSION(2,2), INTENT (IN), OPTIONAL :: src
  INTEGER (kind=long), DIMENSION(2,2), INTENT (IN), OPTIONAL :: dest

  ! filename: is an OPTIONAL STRING OF CHARACTERS of INTENT(IN).
  ! It holds the name of the file to be used
  CHARACTER (len=*), INTENT (IN), OPTIONAL :: filename


  ! file_size: is an OPTIONAL INTEGER of INTENT(IN). It holds the length
  ! of the files to be used if necessary
  INTEGER (kind=long), INTENT (IN), OPTIONAL :: file_size

  ! mode: is an OPTIONAL INTEGER of INTENT(IN). If mode==0, then the
  ! subroutine will firstly try to use a temporary array. If mode==1,
  ! then the subroutine will immediately use temporary files
  INTEGER, INTENT (IN), OPTIONAL :: mode


  ! --------------------------------------
  ! Local variables
  ! --------------------------------------

  ! wtemp: is a REAL allocatable array
  INTEGER (kind=myinteger), DIMENSION (:,:), ALLOCATABLE :: wtemp

  ! lengths to use for wtemp
  INTEGER :: lwtemp1, lwtemp2

  ! Fortran stat parameter
  INTEGER :: stat

  INTEGER :: number_files, mode_copy
  INTEGER (kind=long) :: file_size_copy
  INTEGER, DIMENSION (:), ALLOCATABLE :: units

  integer(long), dimension(2,2) :: src_copy, dest_copy

  info%flag = 0

  ! --------------------------------------
  ! Check input for errors
  ! --------------------------------------

  ! Check w is allocated
  if(.not.allocated(w)) then
    info%flag = zb01_err_w_unalloc
    return
  endif

  src_copy(1,1) = 1
  src_copy(2,1) = size_in(1)
  src_copy(1,2) = 1
  src_copy(2,2) = size_in(2)
  if(present(src)) then
    src_copy(:,:) = src(:,:)
  elseif(present(dest)) then
    src_copy(:,:) = dest(:,:)
  endif
  dest_copy(:,:) = src_copy(:,:)
  if(present(dest)) dest_copy(:,:) = dest(:,:)

  if(minval(src_copy).lt.0 .or. minval(dest_copy).lt.0) &
    info%flag = zb01_err_lkeep_l

  if(maxval(src_copy(:,1)).gt.size_in(1) .or. &
      maxval(src_copy(:,2)).gt.size_in(2)) then
    if(info%flag.eq.zb01_err_lkeep_l) then
      info%flag = zb01_err_lkeep_both
    else
      info%flag = zb01_err_lkeep
    endif
  endif
  if(info%flag.lt.0) return

  ! Note: be careful to only return one error for each dimension
  ! (combination of flags possible if dimensions are differently wrong)
  ! this is to match v1.0.0 behaviour
  if(minval(size_out).lt.1) info%flag = zb01_err_lw_1
  if((maxval(dest_copy(:,1)).gt.size_out(1) .and. size_out(1).ge.1) .or. &
     (maxval(dest_copy(:,2)).gt.size_out(2) .and. size_out(2).ge.1)) then
    if(info%flag.eq.zb01_err_lw_1) then
      info%flag = zb01_err_lw_both
    else
      info%flag = zb01_err_lw
    endif
  endif
  if(info%flag.lt.0) return
  
  if(any(src_copy(2,:)-src_copy(1,:) .ne. &
      dest_copy(2,:)-dest_copy(1,:))) then
    info%flag = zb01_err_src_dest
    return
  endif

  IF (present(filename)) THEN
    IF (len(filename)>400) THEN
      ! Error: len(filename) > 400
      info%flag = zb01_err_filename
      RETURN
    END IF
  END IF

  file_size_copy = 2**22_long
  if (present(file_size)) file_size_copy = file_size
  IF (file_size_copy<2**12) THEN
    ! Error: filesize < 2**12
    info%flag = zb01_err_file_size
    RETURN
  END IF

  mode_copy = 0
  if(present(mode)) mode_copy = mode
  IF (mode_copy<0 .OR. mode_copy>1) THEN
    ! Error: mode out of range
    info%flag = zb01_err_mode
    RETURN
  END IF

  IF (any(src_copy(2,:)-src_copy(1,:)+1.le.0)) THEN

    ! --------------------------------------
    ! No entries need to be copied
    ! --------------------------------------

    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out(1)==size_in(1) .AND. size_out(2)==size_in(2)) RETURN

    ! --------------------------------------
    ! Deallocate w
    ! --------------------------------------
    DEALLOCATE (w,STAT=info%stat)
    IF (info%stat>0) THEN
      info%flag = zb01_err_memory_dealloc
      RETURN
    END IF

    ! --------------------------------------
    ! Reallocate w
    ! --------------------------------------
    ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
    IF (info%stat>0) THEN
      ! --------------------------------------
      ! Allocation of w failed
      ! --------------------------------------
      info%flag = zb01_err_memory_alloc
      RETURN
    END IF

  ELSE

    ! --------------------------------------
    ! Entries need to be copied
    ! --------------------------------------


    ! --------------------------------------
    ! Check whether expansion is required
    ! --------------------------------------
    IF (size_out(1)==size_in(1) .AND. size_out(2)==size_in(2) .and. &
         all(src_copy(:,:).eq.dest_copy(:,:))) RETURN

    ! --------------------------------------
    ! Attempt to allocate temporary array
    ! --------------------------------------
    stat = 0
    ! Size of temporary array
    lwtemp1 = src_copy(2,1) - src_copy(1,1) + 1
    lwtemp2 = src_copy(2,2) - src_copy(1,2) + 1

    ! Allocate temporary array
    IF (mode_copy==0) ALLOCATE (wtemp(lwtemp1,lwtemp2),STAT=stat)
    IF (stat==0 .AND. mode_copy==0) THEN
      ! --------------------------------------
      ! Allocation of temporary array successful
      ! --------------------------------------

      ! --------------------------------------
      ! Copy entries into wtemp
      ! --------------------------------------
      !print *, "copy to wtemp"
      wtemp(1:lwtemp1, 1:lwtemp2) = &
        w(src_copy(1,1):src_copy(2,1), src_copy(1,2):src_copy(2,2))

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
      IF (info%stat>0) THEN
        ! --------------------------------------
        ! Allocation of w not successful so copy temp into files
        ! --------------------------------------
        !print *, "copy wtemp to file"
        call write_to_file(wtemp, lwtemp1+0_long, src_copy, &
          units, number_files, file_size_copy, info, filename)
        if(info%flag.lt.0) return


        ! --------------------------------------
        ! Deallocate wtemp
        ! --------------------------------------
        DEALLOCATE (wtemp,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = zb01_err_memory_dealloc
          RETURN
        END IF

        ! --------------------------------------
        ! Reallocate w
        ! --------------------------------------
        ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
        IF (info%stat>0) THEN
          ! --------------------------------------
          ! Allocation of w to desired size failed
          ! Try size=lwtemp
          ! --------------------------------------

          ALLOCATE (w(lwtemp1,lwtemp2),STAT=info%stat)
          IF (info%stat==0) THEN
            ! --------------------------------------
            ! Reassign size_out
            ! --------------------------------------
            size_out(1) = lwtemp1
            size_out(2) = lwtemp2
            info%flag = zb01_warn_lw

          ELSE
            ! Delete files
            call delete_files(units, number_files, info, filename)
            if(info%flag.lt.0) return

            info%flag = zb01_err_memory_alloc
            RETURN
          END IF
        END IF

        !print *, "restore from file1"
        call read_from_file(w, size_in(1), dest_copy, &
          units, number_files, file_size_copy, info, filename)

        RETURN

      END IF

      ! --------------------------------------
      ! Copy entries
      ! --------------------------------------
      w(dest_copy(1,1):dest_copy(2,1), dest_copy(1,2):dest_copy(2,2)) =&
        wtemp(1:lwtemp1, 1:lwtemp2)

      ! --------------------------------------
      ! Deallocate wtemp
      ! --------------------------------------
      DEALLOCATE (wtemp,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      RETURN

    ELSE
      ! --------------------------------------
      ! Allocation of temporary array not successful or mode = 1
      ! --------------------------------------
      call write_to_file(w, size_in(1), src_copy, &
        units, number_files, file_size_copy, info, filename)
      if(info%flag.lt.0) return

      ! --------------------------------------
      ! Deallocate w
      ! --------------------------------------
      DEALLOCATE (w,STAT=info%stat)
      IF (info%stat>0) THEN
        info%flag = zb01_err_memory_dealloc
        RETURN
      END IF

      ! --------------------------------------
      ! Reallocate w
      ! --------------------------------------
      ALLOCATE (w(size_out(1),size_out(2)),STAT=info%stat)
      IF (info%stat>0) THEN
        ! --------------------------------------
        ! Allocation of w to desired size failed
        ! Try size=lwtemp
        ! --------------------------------------

        ALLOCATE (w(lwtemp1,lwtemp2),STAT=info%stat)
        IF (info%stat==0) THEN
          ! --------------------------------------
          ! Reassign lw
          ! --------------------------------------
          size_out(1) = lwtemp1
          size_out(2) = lwtemp2
          info%flag = zb01_warn_lw

        ELSE
          ! --------------------------------------
          ! Allocation of w fialed - Delete files
          ! --------------------------------------
          call delete_files(units, number_files, info, filename)
          if(info%flag.lt.0) return

          info%flag = zb01_err_memory_alloc
          RETURN
        END IF
      END IF

      call read_from_file(w, size_out(1), dest_copy, &
        units, number_files, file_size_copy, info, filename)


    END IF

  END IF

END SUBROUTINE zb01_resize2_integer

subroutine write_to_file(w, ldw, src, units, number_files, &
      file_size_copy, info, filename)
  integer(long), intent(in) :: ldw ! leading edge of array w
  integer(myinteger), dimension(ldw, *), intent(in) :: w
  integer(long), dimension(2,2), intent(in) :: src
  integer, dimension(:), allocatable, intent(out) :: units
  integer, intent(out) :: number_files
  integer(long), intent(in) :: file_size_copy
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  integer :: iolen
  CHARACTER (402) :: filename_no
  integer(long) :: ncol
  integer(long) :: nrow
  integer(long) :: total_size
  integer :: i, u
  CHARACTER (10) :: ci        ! Filename extension
  logical :: ex
  integer(long) :: idx1       ! first index into w
  integer(long) :: idx2       ! second index into w
  integer(long) :: written    ! number of bytes written to current file
  integer(long) :: len        ! length to write to file

  !print *, "write to file ", src(:,1), ",", src(:,2)

  ! --------------------------------------
  ! Work out number of files required
  ! --------------------------------------
  INQUIRE (iolength=iolen) w(1,1)
  nrow = src(2,1) - src(1,1)
  ncol = src(2,2) - src(1,2)
  total_size = nrow*ncol*iolen

  number_files = (total_size-1)/file_size_copy + 1
  info%files_used = number_files

  ! --------------------------------------
  ! Check files do not exist
  ! --------------------------------------
  IF (present(filename)) THEN
    DO i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      INQUIRE (file=trim(filename_no),exist=ex, iostat=info%iostat)
      IF (info%iostat/=0) THEN
          info%flag = zb01_err_inquire
        return
      ELSE IF (ex) THEN
        info%flag = zb01_err_filename_exists
        return
      END IF
    END DO
  END IF


  ! --------------------------------------
  ! Allocate array for holding unit numbers
  ! --------------------------------------
  i = number_files
  IF (present(filename)) i = 1
  ALLOCATE (units(i),STAT=info%stat)
  IF (info%stat>0) THEN
    info%flag = zb01_err_memory_alloc
    return
  END IF

  ! --------------------------------------
  ! Find unit numbers
  ! --------------------------------------
  call find_units(units, info%iostat)
  if(info%iostat.ne.0) return

  ! --------------------------------------
  ! Open temporary files
  ! --------------------------------------

  IF (present(filename)) THEN
    DO i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no),iostat=info%iostat, &
        err=70,status='new',recl=file_size_copy,form='unformatted', &
        action='readwrite')
      CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
    END DO
  ELSE
    DO i = 1, number_files
      ! write(6,*) 'unit=',units(i)
      u = units(i)
      OPEN (unit=u,iostat=info%iostat,err=70,status='scratch', &
        recl=file_size_copy,form='unformatted',action='readwrite')
    END DO
  END IF

  ! --------------------------------------
  ! Copy entries
  ! --------------------------------------

  idx1 = src(1,1)
  idx2 = src(1,2)
  written = 0
  files: do i = 1, number_files
    if(present(filename)) then
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old',recl=file_size_copy, &
        form='unformatted',action='readwrite', &
        position='rewind')
    else
      u = units(i)
    endif
    do while(file_size_copy - written .gt. 0)
       len = src(2,1) - idx1 + 1 ! amount to reach end of current column
       len = min(len, (file_size_copy-written)/iolen) ! limit to file size
       WRITE (unit=u,iostat=info%iostat,err=90) w(idx1:idx1+len-1,idx2)
       written = written + len*iolen
       idx1 = idx1 + len
       if(idx1.gt.src(2,1)) then
          idx1 = src(1,1)
          idx2 = idx2 + 1
          if(idx2.gt.src(2,2)) then
            if(present(filename)) &
              CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
            exit files
          endif
       endif
    end do
    if(present(filename)) &
      CLOSE (unit=u,iostat=info%iostat,err=80,status='keep')
  end do files

  return

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

  90 info%flag = zb01_err_write
  RETURN
end subroutine write_to_file

!
! This subroutine fills the array units with available unit numbers
! on which files can be opened
!
subroutine find_units(units, iost)
   integer, dimension(:), intent(out) :: units
   integer, intent(inout) :: iost

   integer :: i ! element of units we need to find
   integer :: u ! current unit to try
   logical :: ex ! .true. if unit exists
   logical :: open ! .true. if unit is not open

   iost = 0 ! initialise in case we do no loop iterations

   u = 8 ! unit to start with
   DO i = 1, size(units)
     DO u = u, huge(0)
       IF (u==100 .OR. u==101 .OR. u==102) CYCLE
       INQUIRE (unit=u,iostat=iost,exist=ex, opened=open)
       if(iost.ne.0) return
       IF (ex .AND. .NOT. open) THEN
         units(i) = u
         EXIT
       END IF
     END DO
     u = units(i) + 1
   END DO
end subroutine find_units

subroutine read_from_file(w, ldw, dest, units, number_files, &
      file_size_copy, info, filename)
  integer(long), intent(in) :: ldw
  integer(myinteger), dimension(ldw,*), intent(out) :: w
  integer(long), dimension(2,2), intent(in) :: dest
  integer, dimension(:), intent(in) :: units
  integer, intent(in) :: number_files
  integer(long), intent(in) :: file_size_copy
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  integer :: iolen            ! io length of a single element of w
  CHARACTER (402) :: filename_no
  character (10) :: ci
  integer :: i                ! current file
  integer :: u                ! current unit
  integer(long) :: idx1       ! first index into w
  integer(long) :: idx2       ! second index into w
  integer(long) :: done       ! number of bytes read from current file
  integer(long) :: len        ! length to read from file

  !print *, "read from file ", dest(:,1), ",", dest(:,2)

  INQUIRE (iolength=iolen) w(1,1)

  idx1 = dest(1,1)
  idx2 = dest(1,2)
  done = 0
  files: do i = 1, number_files
    if(present(filename)) then
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old',recl=file_size_copy, &
        form='unformatted',action='readwrite', &
        position='rewind')
    else
      u = units(i)
      rewind(u)
    endif
    do while(file_size_copy - done .gt. 0)
       len = dest(2,1) - idx1 + 1 ! amount to reach end of current column
       len = min(len, (file_size_copy-done)/iolen) ! limit to file size
       READ (unit=u,iostat=info%iostat,err=90) w(idx1:idx1+len-1,idx2)
       done = done + len*iolen
       idx1 = idx1 + len
       if(idx1.gt.dest(2,1)) then
          idx1 = dest(1,1)
          idx2 = idx2 + 1
          if(idx2.gt.dest(2,2)) then
            if(present(filename)) &
              CLOSE (unit=u,iostat=info%iostat,err=80,status='delete')
            exit files
          endif
       endif
    end do
    CLOSE (unit=u,iostat=info%iostat,err=80,status='delete')
  end do files

  return

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

  90 info%flag = zb01_err_read
  RETURN
end subroutine read_from_file

subroutine delete_files(units, number_files, info, filename)
  integer, dimension(:), intent(in) :: units
  integer, intent(in) :: number_files
  type(ZB01_info), intent(inout) :: info
  character(len=*), optional, intent(in) :: filename

  character(402) :: filename_no
  character(10) :: ci
  integer :: i, u

  IF (present(filename)) THEN
    do i = 1, number_files
      IF (number_files/=1) THEN
        WRITE (ci,'(i5)') i
        filename_no = trim(filename) // adjustl(ci)
      ELSE
        filename_no = trim(filename)
      END IF
      u = units(1)
      OPEN (unit=u,file=trim(filename_no), &
        iostat=info%iostat,err=70,status='old', &
        form='unformatted',action='readwrite')
      CLOSE (unit=u,iostat=info%iostat,err=80, &
        status='delete')
    END DO
  ELSE
    DO i = 1, number_files
      u = units(i)
      CLOSE (unit=u,iostat=info%iostat,err=80, &
        status='delete')
    END DO
  END IF

  RETURN

  70 info%flag = zb01_err_open
  RETURN

  80 info%flag = zb01_err_close
  RETURN

end subroutine delete_files

END MODULE hsl_zb01_integer
! COPYRIGHT (c) 2010 Science and Technology Facilities Council
! Original date 18 January 2011. Version 1.0.0
!
! Written by: Jonathan Hogg, John Reid, and Sue Thorne
!
! Version 1.3.0
! For version history, see ChangeLog
!
module hsl_mc69_double
   implicit none

   private
   public :: HSL_MATRIX_UNDEFINED,                             &
      HSL_MATRIX_REAL_RECT, HSL_MATRIX_CPLX_RECT,              &
      HSL_MATRIX_REAL_UNSYM, HSL_MATRIX_CPLX_UNSYM,            &
      HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_CPLX_HERM_PSDEF,   &
      HSL_MATRIX_REAL_SYM_INDEF, HSL_MATRIX_CPLX_HERM_INDEF,   &
      HSL_MATRIX_CPLX_SYM,                                     &
      HSL_MATRIX_REAL_SKEW, HSL_MATRIX_CPLX_SKEW
   public :: mc69_cscl_clean, mc69_verify, mc69_print, mc69_csclu_convert, &
      mc69_coord_convert, mc69_set_values, mc69_csrlu_convert, &
      mc69_cscl_convert, mc69_cscu_convert, mc69_csru_convert, &
      mc69_csrl_convert

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: zero = 0.0_wp

   ! matrix types : real
   integer, parameter :: HSL_MATRIX_UNDEFINED      =  0 ! undefined/unknown
   integer, parameter :: HSL_MATRIX_REAL_RECT      =  1 ! real rectangular
   integer, parameter :: HSL_MATRIX_REAL_UNSYM     =  2 ! real unsymmetric
   integer, parameter :: HSL_MATRIX_REAL_SYM_PSDEF =  3 ! real symmetric pos def
   integer, parameter :: HSL_MATRIX_REAL_SYM_INDEF =  4 ! real symmetric indef
   integer, parameter :: HSL_MATRIX_REAL_SKEW      =  6 ! real skew symmetric

   ! matrix types : complex
   integer, parameter :: HSL_MATRIX_CPLX_RECT      = -1 ! complex rectangular
   integer, parameter :: HSL_MATRIX_CPLX_UNSYM     = -2 ! complex unsymmetric
   integer, parameter :: HSL_MATRIX_CPLX_HERM_PSDEF= -3 ! hermitian pos def
   integer, parameter :: HSL_MATRIX_CPLX_HERM_INDEF= -4 ! hermitian indef
   integer, parameter :: HSL_MATRIX_CPLX_SYM       = -5 ! complex symmetric
   integer, parameter :: HSL_MATRIX_CPLX_SKEW      = -6 ! complex skew symmetric

   ! Error flags
   integer, parameter :: MC69_SUCCESS                =  0
   integer, parameter :: MC69_ERROR_ALLOCATION       = -1
   integer, parameter :: MC69_ERROR_MATRIX_TYPE      = -2
   integer, parameter :: MC69_ERROR_N_OOR            = -3
   integer, parameter :: MC69_ERROR_M_NE_N           = -4
   integer, parameter :: MC69_ERROR_PTR_1            = -5
   integer, parameter :: MC69_ERROR_PTR_MONO         = -6
   integer, parameter :: MC69_ERROR_ROW_BAD_ORDER    = -7
   integer, parameter :: MC69_ERROR_ROW_OOR          = -8
   integer, parameter :: MC69_ERROR_ROW_DUP          = -9
   integer, parameter :: MC69_ERROR_ALL_OOR          = -10
   integer, parameter :: MC69_ERROR_MISSING_DIAGONAL = -11
   integer, parameter :: MC69_ERROR_IMAG_DIAGONAL    = -12
   integer, parameter :: MC69_ERROR_MISMATCH_LWRUPR  = -13
   integer, parameter :: MC69_ERROR_UPR_ENTRY        = -14
   integer, parameter :: MC69_ERROR_VAL_MISS         = -15
   integer, parameter :: MC69_ERROR_LMAP_MISS        = -16


   ! warning flags
   integer, parameter :: MC69_WARNING_IDX_OOR          = 1
   integer, parameter :: MC69_WARNING_DUP_IDX          = 2
   integer, parameter :: MC69_WARNING_DUP_AND_OOR      = 3
   integer, parameter :: MC69_WARNING_MISSING_DIAGONAL = 4
   integer, parameter :: MC69_WARNING_MISS_DIAG_OORDUP = 5

!            Possible error returns:

! MC69_ERROR_ALLOCATION         Allocation error
! MC69_ERROR_MATRIX_TYPE        Problem with matrix_type
! MC69_ERROR_N_OOR              n < 0 or m < 0
! MC69_ERROR_PTR_1              ptr(1) < 1
! MC69_ERROR_PTR_MONO           Error in ptr (not monotonic)
! MC69_ERROR_VAL_MISS           Only one of val_in and val_out is present
! MC69_ERROR_ALL_OOR            All the variables in a column are out-of-range
! MC69_ERROR_IMAG_DIAGONAL      Hermitian case and diagonal not real
! MC69_ERROR_ROW_BAD_ORDER      Entries within a column are not sorted by
!                               increasing row index 
! MC69_ERROR_MISMATCH_LWRUPR    Symmetric, skew symmetric or Hermitian: 
!                               entries in upper and lower
!                               triangles do not match
! MC69_ERROR_MISSING_DIAGONAL   Pos def and diagonal entries missing 
! MC69_ERROR_ROW_OOR            Row contains out-of-range entries      
! MC69_ERROR_ROW_DUP            Row contains duplicate entries         
! MC69_ERROR_M_NE_N             Square matrix and m .ne. n        

!           Possible warnings:

! MC69_WARNING_IDX_OOR          Out-of-range variable indices
! MC69_WARNING_DUP_IDX          Duplicated variable indices
! MC69_WARNING_DUP_AND_OOR      out of range and duplicated indices
! MC69_WARNING_MISSING_DIAGONAL Indefinite case and diagonal entries missing
! MC69_WARNING_MISS_DIAG_OORDUP As MC69_WARNING_MISSING_DIAGONAL, and 
!                               out-of-range and/or duplicates


!!!!!!!!!!!!!!!!!!!!!!!!
! Internal types
!!!!!!!!!!!!!!!!!!!!!!!!
type dup_list
   integer :: src
   integer :: dest
   type(dup_list), pointer :: next => null()
end type dup_list

interface mc69_verify
   module procedure mc69_verify_double
end interface mc69_verify
interface mc69_print
   module procedure mc69_print_double
end interface
interface mc69_cscl_clean
   module procedure mc69_cscl_clean_double
end interface mc69_cscl_clean
interface mc69_set_values
   module procedure mc69_set_values_double
end interface mc69_set_values
interface mc69_cscl_convert
   module procedure mc69_cscl_convert_double
end interface mc69_cscl_convert
interface mc69_cscu_convert
   module procedure mc69_cscu_convert_double
end interface mc69_cscu_convert
interface mc69_csclu_convert
   module procedure mc69_csclu_convert_double
end interface mc69_csclu_convert
interface mc69_csrl_convert
   module procedure mc69_csrl_convert_double
end interface mc69_csrl_convert
interface mc69_csru_convert
   module procedure mc69_csru_convert_double
end interface mc69_csru_convert
interface mc69_csrlu_convert
   module procedure mc69_csrlu_convert_double
end interface mc69_csrlu_convert
interface mc69_coord_convert
   module procedure mc69_coord_convert_double
end interface mc69_coord_convert

contains

!
! To verify that a matrix is in HSL format, or identify why it is not
!
subroutine mc69_verify_double(lp, matrix_type, m, n, ptr, row, flag, more, val)
   integer, intent(in) :: lp ! output unit
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr ! column starts
   integer, dimension(*), intent(in) :: row ! row indices.
     ! Entries within each column must be sorted in order of 
     ! increasing row index. no duplicates and/or out of range entries allowed.
   integer, intent(out) :: flag ! return code
   integer, intent(out) :: more ! futher error information (or set to 0)
   real(wp), dimension(:), optional, intent(in) :: val ! matrix values,if any

   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   logical :: diag ! flag for detection of diagonal
   integer :: i
   integer :: j
   integer :: k
   integer :: last ! last row index
   logical :: lwronly
   integer, dimension(:), allocatable :: ptr2
   integer :: st

   context = 'mc69_verify'
   flag = MC69_SUCCESS

   more = 0 ! ensure more is not undefined.

   ! Check matrix_type
   select case(matrix_type)
   case(0)
      ! Undefined matrix. OK, do nothing
   case(1:4,6)
      ! Real matrix. OK, do nothing
   case(-6:-1)
      ! Complex matrix. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   case default
      ! Out of range value. Issue error
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,lp,flag)
      return
   end select

   ! Check m and n are valid; skip remaining tests if n=0
   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,lp,flag)
      return
   end if
   if(abs(matrix_type).ne.HSL_MATRIX_REAL_RECT .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,lp,flag)
      return
   endif
   if(n == 0 .or. m == 0) return

   ! Check ptr(1) is valid
   if(ptr(1) < 1) then
      more = ptr(1)
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,lp,flag)
      return
   end if

   ! Check ptr is monotonically increasing
   do i = 1, n
      if(ptr(i+1).ge.ptr(i)) cycle ! ptr is monotonically increasing, good.
      flag = MC69_ERROR_PTR_MONO
      more = i + 1
      call mc69_print_flag(context,lp,flag)
      return
   end do

   ! Count number of entries in each row. Also:
   ! * Check ordering of entries
   ! * Check for out-of-range entries
   ! * Check for duplicate entries
   ! * Lack of diagonal entries in real pos. def. case.

   ! ptr2(k+2) is set to hold the number of entries in row k
   allocate(ptr2(m+2), stat=st)
   if(st.ne.0) goto 100
   ptr2(:) = 0

   lwronly = (abs(matrix_type).ne.HSL_MATRIX_UNDEFINED) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_RECT) .and. &
             (abs(matrix_type).ne.HSL_MATRIX_REAL_UNSYM)
   do col = 1, n
      last = -1
      diag = .false.
      do j = ptr(col), ptr(col+1)-1
         k = row(j)
         ! Check out-of-range
         if(k.lt.1 .or. k.gt.m) then
            flag = MC69_ERROR_ROW_OOR
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(lwronly .and. k.lt.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. k.eq.col) then
            flag = MC69_ERROR_UPR_ENTRY
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check duplicate
         if(k.eq.last) then
            flag = MC69_ERROR_ROW_DUP
            more = j-1
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check order
         if(k.lt.last) then
            flag = MC69_ERROR_ROW_BAD_ORDER
            more = j
            call mc69_print_flag(context,lp,flag)
            return
         endif
         ! Check for diagonal
         diag = diag .or. (k.eq.col)
         ! Increase count for row k
         ptr2(k+2) = ptr2(k+2) + 1
         ! Store value for next loop
         last = k
      end do
      ! If marked as positive definite, check if diagonal was present
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF .and. .not.diag) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         more = col
         call mc69_print_flag(context,lp,flag)
         return
      endif
   end do

   if(present(val)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr(j)
            ! Note: column cannot be empty as previously checked pattern ok
            if(real(val(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               more = k
               call mc69_print_flag(context,lp,flag)
               return
            end if 
         end do
      end select
   endif

   return

   100 continue ! Allocation error
   flag = MC69_ERROR_ALLOCATION
   more = st
   call mc69_print_flag(context,lp,flag)
   return

end subroutine mc69_verify_double

!****************************************

!
! Pretty prints a matrix as best it can
!
subroutine mc69_print_double(lp, lines, matrix_type, m, n, ptr, row, val, cbase)
   integer, intent(in) :: lp ! unit to print on
   integer, intent(in) :: lines ! max number of lines to use (ignored if -ive)
   integer, intent(in) :: matrix_type ! type of matrix
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of cols in matrix
   integer, dimension(n+1), intent(in) :: ptr ! column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! row indices
   real(wp), dimension(ptr(n+1)-1), optional, intent(in) :: val ! matrix vals
   logical, optional, intent(in) :: cbase ! is true, rebase for C

   integer :: col, j, k
   integer :: llines
   integer, dimension(:,:), allocatable :: dmat
   character(len=5) :: mfrmt, nfrmt, nefrmt
   character(len=12) :: negfrmt, valfrmt, emptyfrmt
   integer ::  rebase

   if(lp.lt.0) return ! invalid unit

   ! Check if we need to rebase for C consistent output
   rebase = 0
   if(present(cbase)) then
      if(cbase) rebase = 1
   endif

   ! Calculate number of lines to play with
   llines = huge(llines)
   if(lines.gt.0) llines = lines

   ! Print a summary statement about the matrix
   mfrmt = digit_format(m)
   nfrmt = digit_format(n)
   nefrmt = digit_format(ptr(n+1)-1)

   select case(matrix_type)
   case(HSL_MATRIX_UNDEFINED)
      write(lp, "(a)", advance="no") &
         "Matrix of undefined type, dimension "
   case(HSL_MATRIX_REAL_RECT)
      write(lp, "(a)", advance="no") &
         "Real rectangular matrix, dimension "
   case(HSL_MATRIX_REAL_UNSYM)
      write(lp, "(a)", advance="no") &
         "Real unsymmetric matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric positive definite matrix, dimension "
   case(HSL_MATRIX_REAL_SYM_INDEF)
      write(lp, "(a)", advance="no") &
         "Real symmetric indefinite matrix, dimension "
   case(HSL_MATRIX_REAL_SKEW)
      write(lp, "(a)", advance="no") &
         "Real skew symmetric matrix, dimension "
   case default
      write(lp, "(a,i5)") &
         "Unrecognised matrix_type = ", matrix_type
      return
   end select
   write(lp, mfrmt, advance="no") m
   write(lp, "(a)", advance="no") "x"
   write(lp, nfrmt, advance="no") n
   write(lp, "(a)", advance="no") " with "
   write(lp, nefrmt, advance="no") ptr(n+1)-1
   write(lp, "(a)") " entries."

   ! immediate return if m = 0 or n = 0
   if (m == 0 .or. n == 0) return

   if(((present(val) .and. n.lt.10) .or. (.not.present(val) .and. n.lt.24)) &
         .and. m+1.le.llines) then
      ! Print the matrix as if dense
      allocate(dmat(m, n))
      dmat(:,:) = 0
      do col = 1, n
         do j = ptr(col), ptr(col+1) - 1
            k = row(j)
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) then
               dmat(col, k) = -j
            endif
            dmat(k, col) = j
         end do
      end do

      select case(n)
      case(:6)
         valfrmt = "(1x,es12.4)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,12x)"
      case(7)
         valfrmt = "(1x,es10.2)"
         negfrmt = valfrmt
         emptyfrmt = "(1x,10x)"
      case(8:)
         valfrmt = "(1x,es8.2)"
         negfrmt = "(1x,es8.1)"
         emptyfrmt = "(1x,8x)"
      end select

      do k = 1, m
         write(lp,mfrmt,advance="no") k-rebase
         write(lp,"(':')",advance="no")
         if(present(val)) then
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,emptyfrmt,advance="no")
               elseif(dmat(k,j).gt.0) then
                  if(val(dmat(k,j)).gt.zero) then
                     write(lp,valfrmt,advance="no") val(dmat(k,j))
                  else
                     write(lp,negfrmt,advance="no") val(dmat(k,j))
                  endif
               else
                  ! in upper triangle
                  select case(matrix_type)
                  case(HSL_MATRIX_REAL_SYM_INDEF,HSL_MATRIX_REAL_SYM_PSDEF)
                     if(val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") val(-dmat(k,j))
                     endif
                  case(HSL_MATRIX_REAL_SKEW)
                     if(-val(-dmat(k,j)).gt.zero) then
                        write(lp,valfrmt,advance="no") -val(-dmat(k,j))
                     else
                        write(lp,negfrmt,advance="no") -val(-dmat(k,j))
                     endif
                  end select
               endif
            end do
         else ! pattern only
            do j = 1, n
               if(dmat(k,j).eq.0) then 
                  ! nothing here
                  write(lp,"(2x)",advance="no")
               else
                  write(lp,"(1x,'x')",advance="no")
               endif
            end do
         end if
         write(lp,"()")
      end do
   else
      ! Output first x entries from each column
      llines = llines - 1 ! Discount first info line
      if(llines.le.2) return
      write(lp, "(a)") "First 4 entries in columns:"
      llines = llines - 1
      do col = 1, llines
         write(lp, "(a)", advance="no") "Col "
         write(lp, nfrmt, advance="no") col-rebase
         write(lp, "(':')", advance="no")
         do j = ptr(col), min(ptr(col+1)-1, ptr(col)+3)
            write(lp, "(2x)", advance="no")
            write(lp, mfrmt, advance="no") row(j)-rebase
            if(present(val)) then
               write(lp, "(' (',es12.4,')')", advance="no") val(j)
            endif
         end do
         write(lp, "()")
      end do
   endif
end subroutine mc69_print_double

!****************************************

character(len=5) function digit_format(x)
   integer, intent(in) :: x

   integer :: ndigit

   ndigit = log10(real(x)) + 1
   if(ndigit<10) then
      write(digit_format,"('(i',i1,')')") ndigit
   else
      write(digit_format,"('(i',i2,')')") ndigit
   endif
end function digit_format

!****************************************

! Convert a matrix held in lower compressed sparse column format to standard
! HSL format in-place. Duplicates are summed and out-of-range entries are
! remove.  For symmetric, skew-symmetric and Hermitian matrices, entries in 
! the upper triangle are removed and discarded.
!
! Note: this routine is in place! cscl_convert is the out of place version.
subroutine mc69_cscl_clean_double(matrix_type, m, n, ptr, row, flag, &
      val, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m, n ! matrix dimensions
   integer, dimension(n+1), intent(inout) :: ptr ! column pointers
   integer, dimension(*), intent(inout) :: row ! row indices
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(inout) :: val ! values
   integer, optional, intent(out) :: lmap ! size of map
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means val_out(i)=val(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag ! number of columns with a diagonal entry
   integer :: idup ! number of duplicates summed
   integer :: ioor ! number of out-of-range entries
   integer :: j
   integer :: k
   integer :: k1
   integer :: k2
   integer :: l
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer, allocatable :: temp(:) ! work array

   context = 'mc69_cscl_clean'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type) > 1 .and. m /= n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if
 
   if(ptr(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate work array

   allocate(temp(m),stat=st)
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,lp,flag)
      return
   endif
   temp(:) = 0

   ! count duplicates and out-of-range indices
   idup = 0; ioor = 0; idiag = 0
   k = 0 ! last diagonal seen
   l = 1
   do col = 1, n
      if(ptr(col+1).lt.ptr(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      if(abs(matrix_type) > 2) l = col
      do i = ptr(col), ptr(col+1)-1
         j = row(i)
         if(j<l .or. j>m) then
            ioor = ioor + 1
            row(i) = m+1
         else
            if(temp(j)==col) idup = idup + 1
            temp(j) = col
         end if
         if(j.eq.col .and. k<col) then
            idiag = idiag + 1
            k = col
         endif
      end do
   end do
   if(present(ndup)) ndup = idup
   if(present(noor)) noor = ioor
   deallocate(temp,stat=st)
   

   k = 0
   l = ptr(n+1) 
   if(present(map)) then
      deallocate(map,stat=st)
      lmap = l - 1 - ioor + idup
      allocate(map(l-1+idup*2),stat=st)
      if(st /= 0) then
         flag = MC69_ERROR_ALLOCATION
         call mc69_print_flag(context,lp,flag)
         return
      endif
      do i = 1,ptr(n+1)-1
         map(i) = i
      end do
      if(present(val)) then
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:), val=val(k1))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
               val(k) = val(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
                  val(k) = val(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  val(k) = val(k) + val(i)
                  l = l + 2
               end if
            end do
         end do
      else 
         do col = 1, n
            k1 = ptr(col)
            ptr(col) = k + 1
            k2 = ptr(col+1)-1
            if(k2-k1<0) cycle
            call sort (row(k1), k2-k1+1, map=map(k1:))
            ! Squeeze out duplicates and out-of-range indices
            if(row(k1) /= m+1) then
               k = k + 1      
               row(k) = row(k1)
               map(k) = map(k1)
            end if
            do i = k1+1,k2
               if(row(i) == m+1) then
                  exit             
               else if(row(i)>row(i-1)) then
                  k = k + 1
                  row(k) = row(i)
                  map(k) = map(i)
               else
                  map(l) = k
                  map(l+1) = map(i)
                  l = l + 2
               end if
            end do
         end do
      end if
      l = ptr(n+1) - 1
      ptr(n+1) = k + 1
      ! Move duplicate pair in map forward
      do i = 1, idup*2
         map(k+i) = map(l+i)
      end do
   else if(present(val)) then
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1, val=val(k1))
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
            val(k) = val(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
               val(k) = val(i)
            else
               val(k) = val(k) + val(i)
               l = l + 2
            end if
         end do
      end do
      ptr(n+1) = k + 1
   else 
      do col = 1, n
         k1 = ptr(col)
         ptr(col) = k + 1
         k2 = ptr(col+1)-1
         if(k2-k1<0) cycle
         call sort (row(k1), k2-k1+1)
         ! Squeeze out duplicates and out-of-range indices
         if(row(k1) /= m+1) then
            k = k + 1      
            row(k) = row(k1)
         end if
         do i = k1+1,k2
            if(row(i) == m+1) then
               exit             
            else if(row(i)>row(i-1)) then
               k = k + 1
               row(k) = row(i)
            end if
         end do
      end do
      ptr(n+1) = k + 1
   end if

   select case(matrix_type)
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      ! Check for positive diagonal entries
      do j = 1,n
         k = ptr(j)
         ! Note: we're OK if the column is empty - row(k) is still a diagonal
         ! entry, however we must be careful that we've not gone past the
         ! end of the matrix
         if(k.ge.ptr(n+1)) exit
         if(row(k)/=j)then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
         if(present(val))then
            if(val(k)<=zero)then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end if 
      end do
   end select

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

end subroutine mc69_cscl_clean_double

!****************************************

!
! Converts CSC with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_cscl_convert_double(matrix_type, m, n, ptr_in, row_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscl_convert_double

!****************************************

!
! Converts CSC (with lower entries only for symmetric, skew-symmetric and
! Hermitian matrices) to HSL standard format.
! Also used for symmetric, skew-symmetric and Hermitian matrices in upper
! CSR format
!
subroutine mc69_cscl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 or 1, differs for csc/csr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter
   integer :: minidx

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   nullify(dup, duphead)

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! ensure output arrays are not allocated

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)

   idup = 0; ioor = 0; idiag = 0

   allocate(row_out(ptr_in(n+1)-1),stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Allocate map for worst case where all bar one are duplicates
      allocate(map(2*ptr_in(n+1)-2),stat=st)
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            map(k) = multiplier*i
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, map=map(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               call cleanup_dup(duphead)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  allocate(dup,stat=st)
                  if(st.ne.0) goto 100
                  dup%next => duphead
                  duphead => dup
                  dup%src = map(i)
                  dup%dest = k-1
                  cycle
               endif
               if(j.eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               map(k) = map(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            call cleanup_dup(duphead)
            return
         endif
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      allocate(val_out(ptr_in(n+1)-1),stat=st)
      if(st.ne.0) goto 100
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         select case(matrix_type)
         case(HSL_MATRIX_REAL_SKEW)
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = multiplier*val_in(i)
               k = k + 1
            end do
         case default
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.minidx .or. j.gt.m) then
                  ! out of range, ignore
                  ioor = ioor + 1
                  cycle
               endif
               row_out(k) = row_in(i)
               val_out(k) = val_in(i)
               k = k + 1
            end do
         end select
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i, &
               val=val_out(ptr_out(col):k-1))
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, sum then drop from pattern
                  idup = idup + 1
                  val_out(i-1) = val_out(i-1) + val_out(i)
                  cycle
               endif
               if(j.eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               val_out(k) = val_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         endif
      end do
      ptr_out(n+1) = k
   else ! pattern only
      k = 1 ! insert location
      do col = 1, n
         ptr_out(col) = k
         if(ptr_in(col+1).lt.ptr_in(col)) then
            flag = MC69_ERROR_PTR_MONO
            call mc69_print_flag(context,nout,flag)
            return
         endif
         minidx = 1
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) minidx = col
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) minidx = col + 1
         ! Loop over column, copy across while dropping any out of range entries
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.minidx .or. j.gt.m) then
               ! out of range, ignore
               ioor = ioor + 1
               cycle
            endif
            row_out(k) = row_in(i)
            k = k + 1
         end do
         ! Sort entries into order
         i = k - ptr_out(col)
         if(i.eq.0 .and. ptr_in(col+1)-ptr_in(col).ne.0) then
            flag = MC69_ERROR_ALL_OOR
            call mc69_print_flag(context,nout,flag)
            return
         endif
         if(i.ne.0) then
            call sort(row_out(ptr_out(col):k-1), i)
            ! Loop over sorted list and drop duplicates
            i = k-1 ! last entry in column
            k = ptr_out(col)+1 ! insert position
            ! Note: we are skipping the first entry as it cannot be a duplicate
            if(row_out(ptr_out(col)).eq.col) then
               idiag = idiag + 1
            elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            endif
            do i = ptr_out(col)+1, i
               if(row_out(i).eq.row_out(i-1)) then
                  ! duplicate, drop
                  idup = idup + 1
                  cycle
               endif
               if(j.eq.col) idiag = idiag + 1
               row_out(k) = row_out(i)
               k = k + 1
            end do
         elseif(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_PSDEF) then
            flag = MC69_ERROR_MISSING_DIAGONAL
            call mc69_print_flag(context,nout,flag)
            return
         end if
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_cscl_convert_main

!****************************************

!
! Converts CSC with only upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_cscu_convert_double(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: row_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_cscu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_cscu_convert_double

!****************************************

!
! Converts CSC with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csclu_convert_double(matrix_type, n, ptr_in, row_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csclu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, -1, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csclu_convert_double

!****************************************

subroutine mc69_csclu_convert_main(context, multiplier, matrix_type, n, &
      ptr_in, row_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! -1 for upr source, +1 for lwr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: row_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: col ! current column
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: nlwr ! number of strictly lower triangular entries
   integer :: npre
   integer :: nupr ! number of strictly upper triangular entries
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0
   nlwr = 0; nupr = 0

   !
   ! First pass, count number of entries in each row of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   do col = 1, n
      if(ptr_in(col+1).lt.ptr_in(col)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = nlwr + nupr + idiag
      if(ptr_in(col+1).eq.ptr_in(col)) cycle ! no entries in column
      do i = ptr_in(col), ptr_in(col+1)-1
         j = row_in(i)
         if(j.lt.1 .or. j.gt.n) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) then
            ioor = ioor + 1
            cycle
         endif
         if(j.gt.col) then
            nlwr = nlwr + 1
            cycle
         endif
         if(j.lt.col) then
            nupr = nupr + 1
         elseif(j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = col
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(nlwr+nupr+idiag.eq.npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Check number in lower triangle = number in upper triangle
   if(nlwr .ne. nupr) then
      flag = MC69_ERROR_MISMATCH_LWRUPR
      call mc69_print_flag(context,nout,flag)
      return
   endif

   ! Determine column starts for transposed matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle ! ignore oor and lwr triangle entries
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(n+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.ge.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do col = 1, n
            do i = ptr_in(col), ptr_in(col+1)-1
               j = row_in(i)
               if(j.lt.1 .or. j.gt.col) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = col
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      do col = 1, n
         do i = ptr_in(col), ptr_in(col+1)-1
            j = row_in(i)
            if(j.lt.1 .or. j.gt.col) cycle
            if(j.eq.col .and. abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = col
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csclu_convert_main

!****************************************

!
! Converts CSR with only lower entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrl_convert_double(matrix_type, m, n, ptr_in, col_in, ptr_out,&
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup)
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrl_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 1 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csrl_convert_main(context, 1, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrl_convert_double

!****************************************

subroutine mc69_csrl_convert_main(context, multiplier, matrix_type, m, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
   character(50), intent(in) :: context  ! Procedure name (used when printing).
   integer, intent(in) :: multiplier ! 1 for lwr triangle source, -1 for upr
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! matrix dimension
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   integer :: row ! current row
   integer :: col ! current col
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k
   integer :: ne    
   integer :: nout ! output unit (set to -1 if nout not present)
   integer :: npre ! number of out-of-range entries prior to this column
   integer :: st ! stat parameter
   integer :: maxv ! maximum value before out of range

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   if(n < 0 .or. m < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(ptr_in(1) < 1) then
      flag = MC69_ERROR_PTR_1
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix.
   ! Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   k = 0 ! last diagonal found
   maxv = n
   do row = 1, m
      if(ptr_in(row+1).lt.ptr_in(row)) then
         flag = MC69_ERROR_PTR_MONO
         call mc69_print_flag(context,nout,flag)
         return
      endif
      npre = ioor
      if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
      do i = ptr_in(row), ptr_in(row+1)-1
         j = col_in(i)
         if(j.lt.1 .or. j.gt.maxv) then
            ioor = ioor + 1
            cycle
         endif
         if(j.eq.row .and. j.ne.k) then ! diagonal entry (not second in column)
            idiag = idiag + 1
            k = row
         endif
         ptr_out(j+1) = ptr_out(j+1) + 1
      end do
      if(ioor-npre.ne.0 .and. ptr_in(row+1)-ptr_in(row).eq.ioor-npre) then
         ! Column contains only out of range entries
         flag = MC69_ERROR_ALL_OOR
         call mc69_print_flag(context,nout,flag)
         return
      endif
   end do

   ! Check for missing diagonals in pos def case
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
     if(idiag < n) then
         flag = MC69_ERROR_MISSING_DIAGONAL
         call mc69_print_flag(context,nout,flag)
         return
      end if
   end if

   ! Determine column starts for matrix such 
   ! that column i starts at ptr_out(i+1)
   ne = 0
   ptr_out(1) = 1
   do i = 1, n
      ne = ne + ptr_out(i+1)
      ptr_out(i+1) = ptr_out(i) + ptr_out(i+1)
   end do
   do i = n,1,-1
      ptr_out(i+1) = ptr_out(i)
   end do

   !
   ! Second pass, drop entries into place for matrix
   !
   allocate(row_out(ne), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      ! Needs to be this big for worst case: all entries are repeat of same
      allocate(map(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle ! ignore oor and upr triangle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
            map(k) = multiplier*i
         end do
      end do
   elseif(present(val_out)) then
      allocate(val_out(2*(ptr_in(m+1)-1)), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do row = 1, m
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.ge.row) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = multiplier*val_in(i)
            end do
         end do
      case default
         maxv= n
         do row = 1, m
            if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
            if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
            do i = ptr_in(row), ptr_in(row+1)-1
               j = col_in(i)
               if(j.lt.1 .or. j.gt.maxv) cycle
               k = ptr_out(j+1)
               ptr_out(j+1) = k + 1
               row_out(k) = row
               val_out(k) = val_in(i)
            end do
         end do
      end select
   else ! neither val_out nor map present
      maxv = n
      do row = 1, m
         if(abs(matrix_type).ge.HSL_MATRIX_REAL_SYM_PSDEF) maxv = row
         if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW) maxv = row-1
         do i = ptr_in(row), ptr_in(row+1)-1
            j = col_in(i)
            if(j.lt.1 .or. j.gt.maxv) cycle
            k = ptr_out(j+1)
            ptr_out(j+1) = k + 1
            row_out(k) = row
         end do
      end do
   endif

   !
   ! Third pass, removal of duplicates (no sort required by construction)
   !
   if(present(map)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         map(k) = map(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   elseif(present(val_out)) then
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         val_out(k) = val_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               val_out(k-1) = val_out(k-1) + val_out(i)
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else ! neither val_out nor map are present
      k = 1 ! insert position
      do col = 1, n
         i = ptr_out(col)
         ptr_out(col) = k
         if(ptr_out(col+1).eq.i) cycle ! no entries
         ! Move first entry of column forward
         row_out(k) = row_out(i)
         k = k + 1
         ! Loop over remaining entries
         do i = i+1, ptr_out(col+1)-1
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               cycle
            endif
            ! Pull entry forwards
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   endif

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         idup = idup + 1
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup

   return

100 if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
    end if

end subroutine mc69_csrl_convert_main

!****************************************

!
! Converts CSR with upper entries only for symmetric, skew-symmetric and
! Hermitian matrices to HSL standard format
!
subroutine mc69_csru_convert_double(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! number of columns
   integer, dimension(*), intent(in) :: ptr_in ! column pointers on input
   integer, dimension(*), intent(in) :: col_in ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if lp not present)

   context = 'mc69_csru_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_cscl_convert_main(context, -1, matrix_type, n, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csru_convert_double

!****************************************

!
! Converts CSR with both lower and upper entries present to CSC with 
! only lower entries present and entries
! within each column ordered by increasing row index (HSL standard format)
!
subroutine mc69_csrlu_convert_double(matrix_type, n, ptr_in, col_in, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: n ! matrix dimension
   integer, dimension(*), intent(in) :: ptr_in ! row pointers on input
   integer, dimension(*), intent(in) :: col_in ! col indices on input.
      ! These may be unordered within each row and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(*), intent(out) :: ptr_out ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives src: map(i) = j means val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed

   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: nout ! output unit (set to -1 if nout not present)

   context = 'mc69_csrlu_convert'

   nout = -1
   if(present(lp)) nout = lp

   ! Note: have to change this test for complex code
   if(matrix_type < 3 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   call mc69_csclu_convert_main(context, 1, matrix_type, n, &
      ptr_in, col_in, ptr_out, row_out, flag, val_in, val_out, lmap, map, &
      lp, noor, ndup) 
end subroutine mc69_csrlu_convert_double

!****************************************

!
! Converts COOR format to CSC with only lower entries present for
! (skew-)symmetric problems. Entries within each column ordered by increasing
! row index (HSL standard format)
!
subroutine mc69_coord_convert_double(matrix_type, m, n, ne, row, col, ptr_out, &
      row_out, flag, val_in, val_out, lmap, map, lp, noor, ndup) 
   integer, intent(in) :: matrix_type ! what sort of symmetry is there?
   integer, intent(in) :: m ! number of rows in matrix
   integer, intent(in) :: n ! number of columns in matrix
   integer, intent(in) :: ne ! number of input nonzero entries
   integer, dimension(:), intent(in) :: row(ne) ! row indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(in) :: col(ne) ! column indices on input.
      ! These may be unordered within each column and may contain
      ! duplicates and/or out-of-range entries
   integer, dimension(:), intent(out) :: ptr_out(n+1) ! col ptr output
   integer, allocatable, dimension(:), intent(out) :: row_out ! row indices out
      ! Duplicates and out-of-range entries are dealt with and
      ! the entries within each column are ordered by increasing row index.
   integer, intent(out) :: flag ! return code
   real(wp), optional, dimension(*), intent(in) :: val_in ! values input
   real(wp), optional, allocatable, dimension(:) :: val_out
      ! values on output
   integer, optional, intent(out) :: lmap
   integer, optional, allocatable, dimension(:) :: map
      ! map(1:size(val_out)) gives source: map(i) = j means
      ! val_out(i)=val_in(j).
      ! map(size(val_out)+1:) gives pairs: map(i:i+1) = (j,k) means
      !     val_out(j) = val_out(j) + val_in(k)
   integer, optional, intent(in) :: lp ! unit for printing output if wanted
   integer, optional, intent(out) :: noor ! number of out-of-range entries
   integer, optional, intent(out) :: ndup ! number of duplicates summed


   ! Local variables
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i
   integer :: idiag
   integer :: idup
   integer :: ioor
   integer :: j
   integer :: k, l1, l2
   integer :: l
   integer :: ne_new
   integer :: nout ! output unit (set to -1 if lp not present)
   integer :: st ! stat parameter

   type(dup_list), pointer :: dup
   type(dup_list), pointer :: duphead

   nullify(dup, duphead)

   context = 'mc69_coord_convert'

   flag = MC69_SUCCESS

   nout = -1
   if(present(lp)) nout = lp

   ! ---------------------------------------------
   ! Check that restrictions are adhered to
   ! ---------------------------------------------

   ! Note: have to change this test for complex code
   if(matrix_type < 0 .or. matrix_type == 5 .or. matrix_type > 6) then
      flag = MC69_ERROR_MATRIX_TYPE
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(m < 0 .or. n < 0) then
      flag = MC69_ERROR_N_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(abs(matrix_type).ge.HSL_MATRIX_REAL_UNSYM .and. m.ne.n) then
      flag = MC69_ERROR_M_NE_N
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(val_in).neqv.present(val_out)) then
      flag = MC69_ERROR_VAL_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   if(present(map).neqv.present(lmap)) then
      flag = MC69_ERROR_LMAP_MISS
      call mc69_print_flag(context,nout,flag)
      return
   end if

   ! ---------------------------------------------

   ! allocate output and work arrays

   deallocate(row_out,stat=st)
   if(present(val_out)) deallocate(val_out,stat=st)
   if(present(map)) deallocate(map,stat=st)


   ! ---------------------------------------------
   ! check for duplicates and/or out-of-range. check diagonal present.
   ! first do the case where values not present

   idup = 0; ioor = 0; idiag = 0

   !
   ! First pass, count number of entries in each col of the matrix
   ! matrix. Count is at an offset of 1 to allow us to play tricks
   ! (ie. ptr_out(i+1) set to hold number of entries in column i
   ! of expanded matrix).
   ! Excludes out of range entries. Includes duplicates.
   !
   ptr_out(1:n+1) = 0
   do l = 1, ne
      i = row(l)
      j = col(l)
      if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) then
         ioor = ioor + 1
         cycle
      endif

      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SKEW .and. i.eq.j) then
         ioor = ioor + 1
         cycle
      endif
   
      select case (abs(matrix_type))
      case (HSL_MATRIX_REAL_SYM_PSDEF:)
         if(i.ge.j) then
            ptr_out(j+1) = ptr_out(j+1) + 1
         else
            ptr_out(i+1) = ptr_out(i+1) + 1
         end if
      case default
          ptr_out(j+1) = ptr_out(j+1) + 1
      end select
   end do


   ! Determine column starts for transposed expanded matrix such 
   ! that column i starts at ptr_out(i)
   ne_new = 0
   ptr_out(1) = 1
   do i = 2, n+1
      ne_new = ne_new + ptr_out(i)
      ptr_out(i) = ptr_out(i) + ptr_out(i-1)
   end do

   ! Check whether all entries out of range
   if(ne.gt.0 .and. ne_new.eq.0) then
      flag = MC69_ERROR_ALL_OOR
      call mc69_print_flag(context,nout,flag)
      return
   end if

   !
   ! Second pass, drop entries into place for conjugate of transposed
   ! expanded matrix
   !
   allocate(row_out(ne_new), stat=st)
   if(st.ne.0) goto 100
   if(present(map)) then
      if(allocated(map)) deallocate(map,stat=st)
      allocate(map(2*ne), stat=st)
      if(st.ne.0) goto 100
      map(:) = 0
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = -l
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               map(k) = l
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               map(k) = l
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            map(k) = l
         end do
      end select
   elseif(present(val_out)) then
      allocate(val_out(ne_new), stat=st)
      if(st.ne.0) goto 100
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = -val_in(l)
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
               val_out(k) = val_in(l)
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
               val_out(k) = val_in(l)
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
            val_out(k) = val_in(l)
         end do
      end select


   else
      ! neither val_out or map present
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SKEW)
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m .or. i.eq.j) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do

      case(HSL_MATRIX_REAL_SYM_PSDEF, HSL_MATRIX_REAL_SYM_INDEF)
          do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            if(i.ge.j) then
               k=ptr_out(j)
               ptr_out(j) = k+1
               row_out(k) = i
            else
               k=ptr_out(i)
               ptr_out(i) = k+1
               row_out(k) = j
            end if
         end do
      case default
         do l = 1, ne
            i = row(l)
            j = col(l)
            if(j.lt.1 .or. j.gt.n .or. i.lt.1 .or. i.gt.m) cycle
            k=ptr_out(j)
            ptr_out(j) = k+1
            row_out(k) = i
         end do
      end select
   endif

   do j=n,2,-1
      ptr_out(j) = ptr_out(j-1)
   end do
   ptr_out(1) = 1

   !
   ! Third pass, in place sort and removal of duplicates
   ! Also checks for diagonal entries in pos. def. case.
   ! Uses a modified insertion sort for best speed on already ordered data
   !
   idup=0
   if(present(map)) then
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.gt.1) call sort ( row_out(l1:l2),l, map=map(l1:l2) )
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         map(k) = map(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
               ! Duplicate
               idup = idup + 1
               allocate(dup,stat=st)
               if(st.ne.0) goto 100
               dup%next => duphead
               duphead => dup
               dup%src = map(i)
               dup%dest = k-1
               cycle
            endif
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            map(k) = map(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
      lmap = k-1
   else if(present(val_out)) then
      ! ADD
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1) call sort( row_out(l1:l2),l,val=val_out(l1:l2) )
         ! ADD
      end do 

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         val_out(k) = val_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              val_out(k-1) = val_out(k-1)+val_out(i)
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            val_out(k) = val_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   else
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         l = l2-l1+1
         if(l.gt.1)call sort ( row_out(l1:l2),l)
         ! ADD
      end do

      ! work through removing duplicates
      k = 1 ! insert position      
      do j = 1, n
         l1 = ptr_out(j)
         l2 = ptr_out(j+1)-1
         ptr_out(j) = k
         ! sort(row_out(l1:l2)) and permute map(l1:l2) accordingly
         l = l2-l1+1
         if(l.eq.0) cycle ! no entries

         ! Move first entry of column forward
         if(row_out(l1).eq.j) idiag = idiag + 1
         row_out(k) = row_out(l1)
         k = k + 1
         ! Loop over remaining entries
         do i = l1+1, l2
            if(row_out(i).eq.row_out(k-1)) then
              idup = idup + 1
              cycle
            end if
            ! Pull entry forwards
            if(row_out(i).eq.j) idiag = idiag + 1
            row_out(k) = row_out(i)
            k = k + 1
         end do
      end do
      ptr_out(n+1) = k
   
   endif


   ! Check for missing diagonals in pos def and indef cases
   ! Note: change this test for complex case
   if(abs(matrix_type) == HSL_MATRIX_REAL_SYM_PSDEF) then
      do j = 1,n
         if(ptr_out(j).lt.ptr_out(n+1)) then
            if(row_out(ptr_out(j)) .ne. j) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if
         end if
      end do
   end if

   if(present(map)) then
      ! Append duplicates to map
      do while(associated(duphead))
         map(lmap+1) = duphead%dest
         map(lmap+2) = duphead%src
         lmap = lmap + 2
         dup => duphead%next
         deallocate(duphead)
         duphead => dup
      end do
      if(present(val_out)) then
         allocate(val_out(ptr_out(n+1)-1), stat=st)
         if(st.ne.0) goto 100
         call mc69_set_values(matrix_type, lmap, map, val_in, ptr_out(n+1)-1, &
            val_out)
      endif
   endif

   if(present(val_out)) then
      select case(matrix_type)
      case(HSL_MATRIX_REAL_SYM_PSDEF)
         ! Check for positive diagonal entries
         do j = 1,n
            k = ptr_out(j)
            ! ps def case - can't reach here unless all entries have diagonal
            if(real(val_out(k))<=zero) then
               flag = MC69_ERROR_MISSING_DIAGONAL
               call mc69_print_flag(context,nout,flag)
               return
            end if 
         end do
      end select
   endif

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
     if(ioor > 0) flag = MC69_WARNING_IDX_OOR
     if(idup > 0) flag = MC69_WARNING_DUP_IDX
     if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
     if(abs(matrix_type) .ne. HSL_MATRIX_REAL_SKEW) then
         if(idiag < n .and. ioor > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n .and. idup > 0) then
            flag = MC69_WARNING_MISS_DIAG_OORDUP
         else if(idiag < n) then
            flag = MC69_WARNING_MISSING_DIAGONAL
         end if
      endif
      call mc69_print_flag(context,nout,flag)
   end if

   ! Check whether a warning needs to be raised
   if(ioor > 0 .or. idup > 0 .or. idiag < n) then
      if(ioor > 0) flag = MC69_WARNING_IDX_OOR
      if(idup > 0) flag = MC69_WARNING_DUP_IDX
      if(idup > 0 .and. ioor > 0) flag = MC69_WARNING_DUP_AND_OOR
      if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. &
            idiag<n .and. ioor>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n .and.&
            idup>0) then
         flag = MC69_WARNING_MISS_DIAG_OORDUP
      else if(abs(matrix_type).eq.HSL_MATRIX_REAL_SYM_INDEF .and. idiag<n) then
         flag = MC69_WARNING_MISSING_DIAGONAL
      end if
      call mc69_print_flag(context,nout,flag)
   end if

   if(present(noor)) noor = ioor
   if(present(ndup)) ndup = idup
   return

   100 continue
   if(st /= 0) then
      flag = MC69_ERROR_ALLOCATION
      call mc69_print_flag(context,nout,flag)
   end if
   return

end subroutine mc69_coord_convert_double

!*************************************************

!
! This subroutine will use map to translate the values of val to val_out
!
subroutine mc69_set_values_double(matrix_type, lmap, map, val, ne, val_out)
   integer, intent(in) :: matrix_type
   integer, intent(in) :: lmap
   integer, dimension(lmap), intent(in) :: map
   real(wp), dimension(*), intent(in) :: val
   integer, intent(in) :: ne
   real(wp), dimension(ne), intent(out) :: val_out

   integer :: i, j, k

   select case(matrix_type)
   case default
      !
      ! Rectangular, Unsymmetric or Symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + val(k)
      end do
   case(HSL_MATRIX_REAL_SKEW)
      !
      ! Skew symmetric Matrix
      !

      ! First set val_out using first part of map
      do i = 1, ne
         j = abs(map(i))
         val_out(i) = sign(1.0,real(map(i)))*val(j)
      end do

      ! Second examine list of duplicates
      do i = ne+1, lmap, 2
         j = abs(map(i))
         k = abs(map(i+1))
         val_out(j) = val_out(j) + sign(1.0,real(map(i+1)))*val(k)
      end do
   end select
end subroutine mc69_set_values_double

!*************************************************

subroutine mc69_print_flag(context,nout,flag)
   integer, intent (in) :: flag, nout
   character (len=*), optional, intent(in) :: context

   if(nout < 0) return
   if(flag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context), &
         '. Error flag = ', flag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context), &
         '. Warning flag = ', flag
   end if

   select case(flag)
   !
   ! Errors
   !
   case(MC69_ERROR_ALLOCATION)
      write (nout,'(a)') ' Allocation error'
   case(MC69_ERROR_MATRIX_TYPE)
       write (nout,'(a)') ' matrix_type has invalid value'
   case(MC69_ERROR_N_OOR)
      write (nout,'(a)') ' m or n is out-of-range'
   case(MC69_ERROR_ALL_OOR)
      write (nout,'(a)') ' All entries in a column out-of-range'
   case(MC69_ERROR_PTR_MONO)
      write (nout,'(a)') ' ptr not monotonic'
   case(MC69_ERROR_PTR_1)
      write (nout,'(a)') ' ptr(1) < 1'
   case(MC69_ERROR_IMAG_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is not real'
   case(MC69_ERROR_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries are not positive'
   case(MC69_ERROR_VAL_MISS)
      write (nout,'(a)') ' Only one of val and val_out is present'
   case(MC69_ERROR_LMAP_MISS)
      write (nout,'(a)') ' Only one of lmap and map is present'
   case(MC69_ERROR_UPR_ENTRY)
      write (nout,'(a)') ' Entry in upper triangle'
   case(MC69_ERROR_M_NE_N)
      write (nout,'(a)') ' m is not equal to n'
   !
   ! Warnings
   !
   case(MC69_WARNING_IDX_OOR)
      write (nout,'(a)') ' out-of-range indices detected'
   case(MC69_WARNING_DUP_IDX)
      write (nout,'(a)') ' duplicate entries detected'
   case(MC69_WARNING_DUP_AND_OOR)
      write (nout,'(a)') &
         ' out-of-range indices detected and duplicate entries detected'
   case(MC69_WARNING_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is missing'
   case(MC69_WARNING_MISS_DIAG_OORDUP)
      write (nout,'(a)') ' one or more diagonal entries is missing and'
      write (nout,'(a)') ' out-of-range and/or duplicate entries detected'
   end select

end subroutine mc69_print_flag

!************************************************************************

!
!   Sort an integer array by heapsort into ascending order.
!
subroutine sort( array, n, map, val )
   integer, intent(in) :: n       ! Size of array to be sorted
   integer, dimension(n), intent(inout) :: array ! Array to be sorted
   integer, dimension(n), optional, intent(inout) :: map
   real(wp), dimension(n), optional, intent(inout) :: val ! Apply same
      ! permutation to val

   integer :: i
   integer :: temp
   real(wp) :: vtemp
   integer :: root

   if(n.le.1) return ! nothing to do

   !
   ! Turn array into a heap with largest element on top (as this will be pushed
   ! on the end of the array in the next phase)
   !
   ! We start at the bottom of the heap (i.e. 3 above) and work our way
   ! upwards ensuring the new "root" of each subtree is in the correct
   ! location
   root = n / 2
   do root = root, 1, -1
      call pushdown(root, n, array, val=val, map=map)
   end do

   !
   ! Repeatedly take the largest value and swap it to the back of the array
   ! then repeat above code to sort the array
   !
   do i = n, 2, -1
      ! swap a(i) and head of heap a(1)
      temp = array(1)
      array(1) = array(i)
      array(i) = temp
      if(present(val)) then
         vtemp = val(1)
         val(1) = val(i)
         val(i) = vtemp
      endif
      if(present(map)) then
         temp = map(1)
         map(1) = map(i)
         map(i) = temp
      endif
      call pushdown(1,i-1, array, val=val, map=map)
   end do
end subroutine sort

!****************************************

! This subroutine will assume everything below head is a heap and will
! push head downwards into the correct location for it
subroutine pushdown(root, last, array, val, map)
   integer, intent(in) :: root
   integer, intent(in) :: last
   integer, dimension(last), intent(inout) :: array
   real(wp), dimension(last), optional, intent(inout) :: val
   integer, dimension(last), optional, intent(inout) :: map

   integer :: insert ! current insert position
   integer :: test ! current position to test
   integer :: root_idx ! value of array(root) at start of iteration
   real(wp) :: root_val ! value of val(root) at start of iteration
   integer :: root_map ! value of map(root) at start of iteration

   ! NB a heap is a (partial) binary tree with the property that given a
   ! parent and a child, array(child)>=array(parent).
   ! If we label as
   !                      1
   !                    /   \
   !                   2     3
   !                  / \   / \
   !                 4   5 6   7
   ! Then node i has nodes 2*i and 2*i+1 as its children

   if(present(val) .and. present(map)) then ! both val and map
      root_idx = array(root)
      root_val = val(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test);
         val(insert) = val(test);
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
      map(insert) = root_map
   elseif(present(val)) then ! val only, not map
      root_idx = array(root)
      root_val = val(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         val(insert) = val(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
      val(insert) = root_val
   elseif(present(map)) then ! map only, not val
      root_idx = array(root)
      root_map = map(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating mapue up
         array(insert) = array(test)
         map(insert) = map(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root mapue into location found
      array(insert) = root_idx
      map(insert) = root_map
   else ! neither map nor val
      root_idx = array(root)
      insert = root
      test = 2*insert
      do while(test.le.last)
         ! First check for largest child branch to descend
         if(test.ne.last) then
            if(array(test+1).gt.array(test)) test = test + 1
         endif
         if(array(test).le.root_idx) exit ! root gets tested here
         ! Otherwise, move on to next level down, percolating value up
         array(insert) = array(test)
         insert = test
         test = 2*insert
      end do
      ! Finally drop root value into location found
      array(insert) = root_idx
   endif

end subroutine pushdown

!****************************************

subroutine cleanup_dup(duphead)
   type(dup_list), pointer :: duphead ! NB: can't have both intent() and pointer

   type(dup_list), pointer :: dup

   do while(associated(duphead))
      dup => duphead%next
      deallocate(duphead)
      duphead => dup
   end do
end subroutine cleanup_dup

end module hsl_mc69_double
! COPYRIGHT (c) 2006 Council for the Central Laboratory
!                    of the Research Councils
! This package may be copied and used in any application, provided no
! changes are made to these or any other lines.
! Original date 21 February 2006. Version 1.0.0.
! 6 March 2007 Version 1.1.0. Argument stat made non-optional

MODULE HSL_ZD11_double

!  ==========================
!  Sparse matrix derived type
!  ==========================

  TYPE, PUBLIC :: ZD11_type
    INTEGER :: m, n, ne
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: id, type
    INTEGER, ALLOCATABLE, DIMENSION(:) :: row, col, ptr
    REAL ( KIND( 1.0D+0 ) ), ALLOCATABLE, DIMENSION(:) :: val
  END TYPE

CONTAINS

   SUBROUTINE ZD11_put(array,string,stat)
     CHARACTER, allocatable :: array(:)
     CHARACTER(*), intent(in) ::  string
     INTEGER, intent(OUT) ::  stat

     INTEGER :: i,l

     l = len_trim(string)
     if (allocated(array)) then
        deallocate(array,stat=stat)
        if (stat/=0) return
     end if
     allocate(array(l),stat=stat)
     if (stat/=0) return
     do i = 1, l
       array(i) = string(i:i)
     end do

   END SUBROUTINE ZD11_put

   FUNCTION ZD11_get(array)
     CHARACTER, intent(in):: array(:)
     CHARACTER(size(array)) ::  ZD11_get
! Give the value of array to string.

     integer :: i
     do i = 1, size(array)
        ZD11_get(i:i) = array(i)
     end do

   END FUNCTION ZD11_get

END MODULE HSL_ZD11_double


! COPYRIGHT (c) 2007-12 Science and Technology Facilities Council
!           and Jacko Koster (Trondheim, Norway)
!
! History: See ChangeLog
!

! To convert from double to single:
! * s/_double/_single/g
! * Change wp
! * Change: MC64AD, MC64DD, MC64ED, MC64FD, MC64RD, MC64QD, mc64ad, mc64wd
!

module hsl_mc64_double
   use hsl_mc34_double
   use hsl_mc69_double
   use hsl_zd11_double
   implicit none

   private
   public :: mc64_control, mc64_info
   public :: mc64_initialize, mc64_matching

   integer, parameter :: wp = kind(0.0d0)

   type mc64_control
!     real(wp) :: relax = 0d0 ! Relaxes matching
      integer :: lp = 6       ! Unit for error messages
      integer :: wp = 6       ! Unit for warning messages
      integer :: sp = -1      ! Unit for statistical output
      integer :: ldiag = 2    ! Controls level of diagnostic output
      integer :: checking = 0 ! Control for checking input
   end type mc64_control

   type mc64_info
      integer :: flag   ! Flags success or failure case
      integer :: more    ! More information on failure
      integer :: strucrank ! Structural rank
      integer :: stat    ! STAT value after allocate failure
   end type mc64_info

   interface mc64_matching
      module procedure mc64_matching_zd11_double
      module procedure mc64_matching_hslstd_double
   end interface mc64_matching

contains

   subroutine mc64_initialize(control)
      type(mc64_control), intent(out) :: control

      type(mc64_control) :: default

!       control%relax = default%relax
        control%lp = default%lp
        control%wp = default%wp
        control%sp = default%sp
        control%ldiag = default%ldiag
        control%checking = default%checking
    end subroutine mc64_initialize

   subroutine mc64_matching_zd11_double(job,matrix,control,info,perm, &
         scale)
      integer, intent(in) :: job ! Control parameter for algorithm choice
      type(zd11_type), intent(in) :: matrix
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(matrix%m + matrix%n)
      real(wp), optional, intent(out) :: scale(matrix%m + matrix%n)

      integer :: matrix_type

      matrix_type = HSL_MATRIX_UNDEFINED

      if (allocated(matrix%id)) then
        if (matrix%id(1) == 'S' .or. matrix%id(1) == 's') then
          matrix_type = HSL_MATRIX_REAL_SYM_INDEF
        endif
      endif

      if (control%checking == 0) then
        if (matrix%n .lt. 1) then
          info%flag = -2
          info%more = matrix%n
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%n out-of-range', info%more
          return
        endif

        if(matrix%ne .ne. matrix%ptr(matrix%n+1)-1) then
          info%flag = -3
          info%more = matrix%ptr(matrix%n+1)-1
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of ptr(n+1)-1!=ne', info%more
          return
        endif
      endif

      call mc64_matching_hslstd_double(job,matrix_type,matrix%m, &
         matrix%n,matrix%ptr,matrix%row,matrix%val,control,info,perm,scale)
   end subroutine mc64_matching_zd11_double

   subroutine mc64_matching_hslstd_double(job,matrix_type,m,n,ptr,row, &
         val,control,info,perm,scale)
      integer, intent(in) :: job ! Control parameter for algorithm choice
      integer, intent(in) :: matrix_type
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(*), intent(in) :: row
      real(wp), dimension(*), intent(in) :: val
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(m + n)
      real(wp), optional, intent(out) :: scale(m + n)

      integer, allocatable :: iw(:)
      real(wp), allocatable :: dw(:)
      real(wp) :: cntl(10)
      integer :: i,j,k,ne,lliw,liw,lldw,ldw,ndiag,num,stat &
                 ,icntl(10),info64(10)
      real(wp), parameter :: zero = 0.0
      logical sym

      stat = 0

! Check input data if checking is requested
      if (control%checking == 0) then
        if (job .lt. 1 .or. job .gt. 5) then
          info%flag = -1
          info%more = job
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of job out-of-range', info%more
          return
        endif

        if (n .lt. 1) then
          info%flag = -2
          info%more = n
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%n out-of-range', info%more
          return
        endif

        if (ptr(n+1)-1 .lt. 1) then
          info%flag = -3
          info%more = ptr(n+1)-1
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of ptr(n+1)-1 out-of-range', info%more
          return
        endif

        if (m .lt. n) then
          info%flag = -4
          info%more = m
          if (control%ldiag>0 .and. control%lp>0 ) &
              write (control%lp,'(/a,i5/a,i5)') &
              'Error return from MC64 with info%flag =', &
              info%flag, &
              'Value of matrix%m less than matrix%n', info%more
          return
        endif

! end of checking
      endif

      ne = ptr(n+1)-1

      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%sp

! Action if matrix is symmetric
      sym = (matrix_type.ge.HSL_MATRIX_REAL_SYM_PSDEF)
      if(sym .and. job.eq.5) then
         call sym_maxprod_match(n,ptr,row,val,control,info,perm,scale)
         return
      endif
      if (sym) then
        ndiag = 0
! Matrix is symmetric ... check only lower triangle supplied
        if (control%checking == 0) then
          do j = 1,n
            do k = ptr(j),ptr(j+1)-1
              if (row(k) .lt. j) then
                info%flag = -8
                info%more = j
                if (control%ldiag>0 .and. control%lp>0 ) &
                  write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                   info%flag, &
       'Input symmetric matrix has entries in upper triangle in column', &
                   info%more
                return
              endif
              if (row(k) == j) ndiag = ndiag + 1
            enddo
          enddo
        else
          do j = 1,n
            do k = ptr(j),ptr(j+1)-1
              if (row(k) == j) ndiag = ndiag + 1
            enddo
          enddo
        endif
! Reset ne for the expanded symmetric matrix
        ne = 2*ne - ndiag
      endif


      if (job == 1) liw = 4*n + m
      if (job == 2) liw = 2*n + 2*m
      if (job == 3) liw = 8*n+2*m+ne
      if (job == 4) liw = 3*n+2*m
      if (job == 5) liw = 3*n+2*m
!     if (job == 6) liw = 3*n+2*m+ne

      if (job == 1) ldw = 1
      if (job == 2) ldw = m
      if (job == 3) ldw = ne
      if (job == 4) ldw = 2*m+ne
      if (job == 5) ldw = 2*m+n+ne
!     if (job == 6) ldw = 3*m+n+ne

! Expand matrix if symmetric
      if (sym) then
        lliw = liw+n+1+2*ne
        lldw = ldw+2*ne
        allocate (iw(lliw),stat=stat)
        if (stat/=0) go to 100
        allocate (dw(lldw),stat=stat)
        if (stat/=0) go to 100

        iw(n+2:n+1+ptr(n+1)-1) = row(1:ptr(n+1)-1)
        dw(1:ptr(n+1)-1) = val(1:ptr(n+1)-1)
        iw(1:n+1) = ptr(1:n+1)
        !call mc34ad(n,iw(n+2),iw,.true.,dw,iw(n+2+2*ne))
        call mc34_expand(n,iw(n+2:),iw,iw(n+2+2*ne:),a=dw)
      else
        allocate (iw(liw),stat=stat)
        if (stat/=0) go to 100
        allocate (dw(ldw),stat=stat)
        if (stat/=0) go to 100
      endif

      icntl(4) = control%checking
      icntl(5:10) = 0

!     if (m==n .and. control%relax == zero .and. job .le. 5) then
      if (m==n) then
! Needed because no ldiag in mc64
        if (control%ldiag .lt. 1) then
          icntl(1) = -1
          icntl(2) = -1
          icntl(3) = -1
        endif
        if (control%ldiag .eq. 1) then
          icntl(2) = -1
          icntl(3) = -1
        endif
        if (control%ldiag .eq. 2) then
          icntl(3) = -1
        endif
        if (sym)  then
! Call HSL F77 code
          call mc64ad(job,n,ne,iw,iw(n+2),dw,num,perm, &
                      liw,iw(n+2+2*ne),ldw,dw(2*ne+1),icntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info64(1) == 2) info%flag = -9
          if (info%flag .lt. 0) go to 80
! Note that in this case m=n
! Set column permutation and row permutation to be the same
          perm(m+1:m+n) = perm(1:m)
        else
          call mc64ad(job,n,ne,ptr,row,val,num,  &
                      perm(m+1),liw,iw,ldw,dw,icntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info%flag .lt. 0) go to 80
! Set row permutation to identity
          do i = 1,m
            perm(i) = i
          enddo
! Invert column permutation
          do i = 1,n
            iw(abs(perm(m+i))) = sign(i,perm(m+i))
          enddo
          do i = 1,n
            perm(m+i) = iw(i)
          enddo
        endif
      else
        icntl(5) = control%ldiag
        icntl(6:10) = 0
!       cntl(1) = control%relax
        cntl(1) = zero
        call mc64_hsl_a(job,m,n,ne,ptr,row,val, &
                  num,perm,liw,iw,ldw,dw,icntl,cntl,info64)
          info%flag   = info64(1)
          info%more   = info64(2)
          if (info%flag .lt. 0) go to 80
! Set column permutation to identity
        do i = 1,n
          perm(m+i) = i
        enddo
      endif

      info%strucrank  = num

      if (present(scale) .and. job == 5) then
      !!! Commented as handled by seperate subroutine now
      !!!  if (sym) then
      !!!    scale(1:m) = (dw(2*ne+1:2*ne+n)+dw(2*ne+n+1:2*ne+2*n))/2
      !!!    scale(m+1:m+n) = scale(1:m)
      !!!  else
          scale(1:m) = dw(1:m)
          scale(m+1:m+n) = dw(m+1:m+n)
      !!!  endif
      endif

  80  deallocate (iw, stat=stat)
      if (stat/=0) go to 100
      deallocate (dw, stat=stat)
      if (stat/=0) go to 100

      return

  100 info%flag = -5
      info%stat = stat
      if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a,i5/a,i5)') &
         'Error return from MC64 with info%flag =', info%flag, &
         'Allocate or deallocate failed with STAT=',stat

   end subroutine mc64_matching_hslstd_double

   subroutine sym_maxprod_match(n,ptr,row,val,control,info,perm,scale)
      integer, intent(in) :: n
      integer, dimension(n+1), intent(in) :: ptr
      integer, dimension(*), intent(in) :: row
      real(wp), dimension(*), intent(in) :: val
      type(mc64_control), intent(in) :: control
      type(mc64_info), intent(out) :: info
      integer, intent(out) :: perm(2*n)
      real(wp), optional, intent(out) :: scale(2*n)

      integer, allocatable :: ptr2(:), row2(:), iw(:), new_to_old(:), &
         old_to_new(:), cperm(:)
      real(wp), allocatable :: val2(:), dw(:), cmax(:), cscale(:)
      real(wp) :: colmax
      integer :: i,j,k,ne,num,stat,nn,j1,j2,jj
      real(wp), parameter :: zero = 0.0

      stat = 0
      ne = ptr(n+1)-1

! Matrix is symmetric ... check only lower triangle supplied
      if (control%checking == 0) then
        do j = 1,n
          do k = ptr(j),ptr(j+1)-1
            if (row(k) .lt. 1 .or. row(k).gt.n) then
              info%flag = -6
              info%more = j
              if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                  info%flag, &
                  'Input symmetric matrix has entries with invald row index &
                  &in position', &
                  info%more
              return
            endif
            if (row(k) .lt. j) then
              info%flag = -8
              info%more = j
              if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
                  'Error return from MC64 with info%flag =', &
                   info%flag, &
                   'Input symmetric matrix has entries in upper triangle in &
                   &position', &
                   info%more
              return
            endif
          enddo
        enddo
      endif
! Reset ne for the expanded symmetric matrix
      ne = 2*ne

! Expand matrix, drop explicit zeroes and take log absolute values
      allocate (ptr2(n+1), row2(ne), val2(ne), &
                iw(5*n), dw(2*n), cmax(n), stat=stat)
      if (stat/=0) then
         info%flag = -5
         info%stat = stat
         if (control%ldiag>0 .and. control%lp>0 ) &
             write (control%lp,'(/a,i5/a,i5)') &
            'Error return from MC64 with info%flag =', info%flag, &
            'Allocate or deallocate failed with STAT=',stat
         return
      endif

      if(control%checking==0) then
         ! checking enabled, look for duplicates
         iw(1:n) = 0
         do i = 1, n
           do j = ptr(i), ptr(i+1)-1
             if(iw(row(j)).ge.i) then
               ! duplicate detected
               info%flag = -7
               info%more = j
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a,i5/a,i5)') &
                   'Error return from MC64 with info%flag =', &
                    info%flag, &
                    'Input symmetric matrix has duplicate entries in &
                    &position', &
                    info%more
               return
             endif
             iw(row(j)) = i
           end do
         end do
      endif

      k = 1
      do i = 1, n
         ptr2(i) = k
         do j = ptr(i), ptr(i+1)-1
            if(val(j).eq.zero) cycle
            row2(k) = row(j)
            val2(k) = abs(val(j))
            k = k + 1
         end do
         ! Following log is seperated from above loop to expose expensive
         ! log operation to vectorization.
         val2(ptr2(i):k-1) = log(val2(ptr2(i):k-1))
      end do
      ptr2(n+1) = k
      call mc34_expand(n, row2, ptr2, iw, a=val2)

! Compute column maximums
      do i = 1, n
         colmax = maxval(val2(ptr2(i):ptr2(i+1)-1))
         cmax(i) = colmax
         val2(ptr2(i):ptr2(i+1)-1) = colmax - val2(ptr2(i):ptr2(i+1)-1)
      end do

      call mc64wd(n,ne,ptr2,row2,val2,perm,num,iw(1),iw(n+1), &
                  iw(2*n+1),iw(3*n+1),iw(4*n+1),dw(1),dw(n+1))

      info%flag   = 0
      info%strucrank  = num

      if(num.eq.n) then ! Full rank
! Note that in this case m=n
! Set column permutation and row permutation to be the same
         perm(n+1:n+n) = perm(1:n)


         if (present(scale)) then
             scale(1:n) = (dw(1:n)+dw(n+1:2*n)-cmax(1:n))/2
             scale(n+1:n+n) = scale(1:n)
         endif
         return
      endif

      ! If we reach this point then structually rank deficient:
      ! Build a full rank submatrix and call mc64wd on it.
      info%flag = 1 ! structually singular warning

      allocate(old_to_new(n),new_to_old(n),cperm(n),stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = -5
         info%stat = stat
         if (control%ldiag>0 .and. control%lp>0 ) &
             write (control%lp,'(/a,i5/a,i5)') &
            'Error return from MC64 with info%flag =', info%flag, &
            'Allocate or deallocate failed with STAT=',stat
         return
      end if

      j = num + 1
      k = 0
      do i = 1,n
         if (perm(i) < 0) then
            ! row i is not part of the matching
            old_to_new(i) = -j
            j = j + 1
         else
            k = k + 1
            ! old_to_new(i) holds the new index for variable i after
            ! removal of singular part and new_to_old(k) is the
            ! original index for k
            old_to_new(i) = k
            new_to_old(k) = i
         end if
      end do

      ! Overwrite ptr2, row2 and val2
      ne = 0
      k = 0
      ptr2(1) = 1
      j2 = 1
      do i = 1,n
         j1 = j2
         j2 = ptr2(i+1)
         ! skip over unmatched entries
         if (perm(i) < 0) cycle
         k = k + 1
         do j = j1,j2-1
            jj = row2(j)
            if (perm(jj) < 0) cycle
            ne = ne + 1
            row2(ne) = old_to_new(jj)
            val2(ne) = val2(j)
         end do
         ptr2(k+1) = ne + 1
      end do
      ! nn is order of non-singular part.
      nn = k
      call mc64wd(nn,ne,ptr2,row2,val2,cperm,num, &
                  iw(1),iw(nn+1),iw(2*nn+1),iw(3*nn+1),iw(4*nn+1), &
                  dw(1),dw(nn+1))

      if(present(scale)) then
          do i = 1,n
             j = old_to_new(i)
             if (j < 0) then
                scale(i) = -huge(scale)
             else
               ! Note: we need to subtract col max using old matrix numbering
               scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
            end if
         end do
      endif

      perm(1:n) = -1
      do i = 1,nn
         j = cperm(i)
         perm(new_to_old(i)) = j
      end do

      do i = 1, n
         if(perm(i).eq.-1) then
            perm(i) = old_to_new(i)
         endif
      end do

      perm(n+1:n+n) = perm(1:n)

      ! Apply Duff and Pralet correction to unmatched row scalings
      if(present(scale)) then
         allocate(cscale(n), stat=stat)
         if(stat/=0) then
            info%flag = -5
            info%stat = stat
            if (control%ldiag>0 .and. control%lp>0 ) &
                write (control%lp,'(/a,i5/a,i5)') &
               'Error return from MC64 with info%flag =', info%flag, &
               'Allocate or deallocate failed with STAT=',stat
            return
         endif
         ! For columns i not in the matched set I, set
         !     s_i = 1 / (max_{k in I} | a_ik s_k |)
         ! with convention that 1/0 = 1
         cscale(1:n) = scale(1:n)
         do i = 1,n
            do j = ptr(i), ptr(i+1)-1
               k = row(j)
               if(cscale(i).eq.-huge(scale).and.cscale(k).ne.-huge(scale)) then
                  ! i not in I, k in I
                  scale(i) = max(scale(i), log(abs(val(j)))+scale(k))
               endif
               if(cscale(k).eq.-huge(scale).and.cscale(i).ne.-huge(scale)) then
                  ! k not in I, i in I
                  scale(k) = max(scale(k), log(abs(val(j)))+scale(i))
               endif
            end do
         end do
         do i = 1,n
            if(cscale(i).ne.-huge(scale)) cycle ! matched part
            if(scale(i).eq.-huge(scale)) then
               scale(i) = zero
            else
               scale(i) = -scale(i)
            endif
         end do
      endif

   end subroutine sym_maxprod_match

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_A (JOB, M, N, NE, IP, IRN, A, NUM, PERM, LIW,&
      IW, LDW, DW, ICNTL, CNTL, INFO)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! Purpose
! =======
!
! This subroutine attempts to find a permutation for an MxN, M>=N,
! sparse matrix A = {a_ij} that makes the permuted matrix have N
! entries on its diagonal.
! If the matrix is structurally nonsingular, the subroutine optionally
! returns a permutation that maximizes the smallest element on the
! diagonal, maximizes the sum of the diagonal entries, or maximizes
! the product of the diagonal entries of the permuted matrix.
! For the latter option, the subroutine also finds scaling factors that
! may be used to scale the matrix so that the nonzero diagonal entries
! of the permuted matrix are one in absolute value and all the
! off-diagonal entries are less than or equal to one in absolute value.
! The natural logarithms of the scaling factors u(i), i=1..M, for the
! rows and v(j), j=1..N, for the columns are returned so that the
! scaled matrix B = {b_ij} has entries b_ij = a_ij * EXP(u_i + v_j).
! The scaling factors are returned by this subroutine, but the actual
! scaling of the matrix has to be performed by the calling program.
!
! Parameters
! ==========
!
      INTEGER NICNTL, NCNTL, NINFO
      PARAMETER (NICNTL = 10, NCNTL = 10, NINFO = 10)

      INTEGER JOB, M, N, NE, NUM, LIW, LDW
      INTEGER IP (N + 1), IRN (NE), PERM (M), IW (LIW)
      INTEGER ICNTL(NICNTL), INFO (NINFO)
      REAL(WP) A (NE)
!
! JOB is an INTEGER variable which must be set by the user to control
! the action. It is not altered by the subroutine.
! Possible values for JOB are:
!   1 Compute a column permutation of the matrix so that the
!     permuted matrix has as many entries on its diagonal as possible.
!     The values on the diagonal are of arbitrary size. HSL subroutine
!     MC21A/MC64Z is used for this. See [1].
!   2 Compute a column permutation of the matrix so that the smallest
!     value on the diagonal of the permuted matrix is maximized.
!     See [3].
!   3 Compute a column permutation of the matrix so that the smallest
!     value on the diagonal of the permuted matrix is maximized.
!     The algorithm differs from the one used for JOB = 2 and may
!     have quite a different performance. See [2].
!   4 Compute a column permutation of the matrix so that the sum
!     of the diagonal entries of the permuted matrix is maximized.
!     See [3].
!   5 Compute a column permutation of the matrix so that the product
!     of the diagonal entries of the permuted matrix is maximized
!     and vectors to scale the matrix so that the nonzero diagonal
!     entries of the permuted matrix are one in absolute value and
!     all the off-diagonal entries are less than or equal to one in
!     absolute value. See [3].
!  Restriction: 1 <= JOB <= 5.
!
! M is an INTEGER variable which must be set by the user to the
!   number of rows of the matrix A. It is not altered by the
!   subroutine. Restriction: M >= N.
!
! N is an INTEGER variable which must be set by the user to the
!   number of columns of the matrix A. It is not altered by the
!   subroutine. Restriction: N >= 1.
!
! NE is an INTEGER variable which must be set by the user to the
!   number of entries in the matrix. It is not altered by the
!   subroutine. Restriction: NE >= 1.
!
! IP is an INTEGER array of length N+1.
!   IP(J), J=1..N, must be set by the user to the position in array IRN
!   of the first row index of an entry in column J. IP(N+1) must be set
!   to NE+1. It is not altered by the subroutine.
!
! IRN is an INTEGER array of length NE.
!   IRN(K), K=1..NE, must be set by the user to hold the row indices of
!   the entries of the matrix. Those belonging to column J must be
!   stored contiguously in the positions IP(J)..IP(J+1)-1. The ordering
!   of the row indices within each column is unimportant. Repeated
!   entries are not allowed. The array IRN is not altered by the
!   subroutine.
!
! A is a REAL array of length NE.
!   The user must set A(K), K=1..NE, to the numerical value of the
!   entry that corresponds to IRN(K).
!   It is not used by the subroutine when JOB = 1.
!   It is not altered by the subroutine.
!
! NUM is an INTEGER variable that need not be set by the user.
!   On successful exit, NUM will be the number of entries on the
!   diagonal of the permuted matrix.
!   If NUM < N, the matrix is structurally singular.
!
! PERM is an INTEGER array of length M that need not be set by the
!   user. On successful exit, PERM can be interpreted in any of the
!   following ways:
!
!   1. If M=N, PERM contains the column permutation.
!      Column PERM(I) of the original matrix is column I in the
!      permuted matrix, I=1..N.
!      (This was the definition of parameter CPERM in versions of
!      MC64AD before version 1.2b)
!
!   2. If M>=N, PERM contains the row permutation.
!      Row I of the original matrix is row ABS(PERM(I)) in the
!      permuted matrix, I=1..M.
!      The rows where PERM(I) is positive constitute an N by N matrix
!      the scaled version of which has ones on the diagonal.
!
! LIW is an INTEGER variable that must be set by the user to
!   the dimension of array IW. It is not altered by the subroutine.
!   Restriction:
!     JOB = 1 :  LIW >=  4N +  M
!     JOB = 2 :  LIW >=  2N + 2M
!     JOB = 3 :  LIW >=  8N + 2M + NE
!     JOB = 4 :  LIW >=  3N + 2M
!     JOB = 5 :  LIW >=  3N + 2M
!
! IW is an INTEGER array of length LIW that is used for workspace.
!
! LDW is an INTEGER variable that must be set by the user to the
!   dimension of array DW. It is not altered by the subroutine.
!   Restriction:
!     JOB = 1 :  LDW not used
!     JOB = 2 :  LDW >=      M
!     JOB = 3 :  LDW >=          NE
!     JOB = 4 :  LDW >=     2M + NE
!     JOB = 5 :  LDW >= N + 2M + NE
!
! DW is a REAL array of length LDW used for workspace.
!   If JOB = 5, on return, DW(i) contains u_i, i=1..M, and
!   DW(M+j) contains v_j, j=1..N.
!
! ICNTL is an INTEGER array of length NICNTL.
!   Its components control the output of MC64AD and must be set by the
!   user before calling MC64AD. They are not altered by the subroutine.
!
!   ICNTL(1) must be set to specify the output stream for
!   error messages. If ICNTL(1) < 0, messages are suppressed.
!
!   ICNTL(2) must be set by the user to specify the output stream for
!   warning messages. If ICNTL(2) < 0, messages are suppressed.
!
!   ICNTL(3) must be set by the user to specify the output stream for
!   diagnostic messages. If ICNTL(3) < 0, messages are suppressed.
!
!   ICNTL(4) must be set by the user to a value other than 0 to avoid
!   checking of the input data.  Setting ICNTL(4) to any
!   other will avoid the checks but is likely to cause problems
!   later if out-of-range indices or duplicates are present.
!   The user should set ICNTL(4) nonzero, if the data is known not
!   to contain such problems. The code will exhibit undefined
!   behaviour in case data checking is not done and the
!   input data does not satisfy the restrictions as listed
!   elsewhere.
!
!   ICNTL(5) must be set by the user to control the printing of
!   diagnostic messages.
!   If ICNTL(5) <= 0, no messages are output.
!   If ICNTL(5) = 1, only error messages are output.
!   If ICNTL(5) = 2, error and warning messages output.
!   If ICNTL(5) = 3, as for 2 plus scalar parameters, the first 
!   ten entries of array parameters, and the control parameters on 
!   the first entry.
!   If ICNTL(5) > 3, full data will be printed on entry and exit.
!
! CNTL is a REAL array of length NCNTL.
!   Its components control the output of MC64AD and must be set by the
!   user before calling MC64AD. They are not altered by the subroutine.
!
!   CNTL(1) must be set to specify the relaxation parameter.
!   It is used by MC64 only if JOB = 3,4,5.
!   It must be set to a non-negative value (usually close to zero).
!   If CNTL(1) < 0.0, it is treated as 0.0.
!
!   CNTL(1) is a relaxation parameter. A positive value will lead to
!   matchings computed by MC64AD that are not optimal/maximal in some
!   sense but only nearly so. However, these non-optimal matchings are
!   often computed more quickly. Appropriate values for CNTL(1) are
!   problem dependent but usually slightly larger than 0.0.
!
!
! INFO is an INTEGER array of length NINFO which need not be set by the
!   user. INFO(1) is set non-negative to indicate success. A negative
!   value is returned if an error occurred, a positive value if a
!   warning occurred. INFO(2) holds further information on the error.
!   On exit from the subroutine, INFO(1) will take one of the
!   following values:
!    0 : successful entry (for structurally nonsingular matrix).
!   +1 : successful entry (for structurally singular matrix).
!   +2 : the returned scaling factors are large and may cause
!        overflow when used to scale the matrix.
!        (For JOB = 4,5 entries only.)
!   +4 : CNTL(1) is negative and treated as zero.
!   -1 : JOB < 1 or JOB > 5.  Value of JOB held in INFO(2).
!   -2 : N < 1.  Value of invalid N held in INFO(2).
!   -3 : NE < 1.  Value of NE held in INFO(2).
!   -4 : M < N. Value of M held in INFO(2).
!   -6 : entries are found whose row indices are out of range. INFO(2)
!        contains the position in arrays A/IRN in which first entry is found.
!        (This value can be returned only if ICNTL(4) was set to zero.)
!   -7 : repeated entries are found. INFO(2) contains the position in arrays
!        A/IRN in which first entry is found.
!        (This value can be returned only if ICNTL(4) was set to zero.)
!
!   A return with one of the values INFO(1)=+3,+5,+6,+7 is also possible
!   These values are combinations of the above warnings (+1,+2,+4) and
!   correspond to the sum of the constituent warnings.
!
!   INFO(3) to INFO(NINFO) are not currently used and are set to zero
!        by the routine.
!
      REAL(WP) DW (LDW), CNTL (NCNTL)
! References:
!  [1] I. S. Duff, (1981),
!      "Algorithm 575. Permutations for a zero-free diagonal",
!      ACM Trans. Math. Software 7(3), 387-390.
!  [2] I. S. Duff and J. Koster, (1998),
!      "The design and use of algorithms for permuting large
!      entries to the diagonal of sparse matrices",
!      SIAM J. Matrix Anal. Appl., vol. 20, no. 4, pp. 889-901.
!  [3] I. S. Duff and J. Koster, (2001),
!      "On algorithms for permuting large entries to the diagonal
!      of sparse matrices",
!      SIAM J. Matrix Anal. Appl., vol. 22, no. 4, pp. 973-996.

! Local variables and parameters
      INTEGER I, J, K, WARN1, WARN2, WARN4
      REAL(WP) FACT, RINF
      REAL(WP), PARAMETER :: ZERO = 0.0_WP
      REAL(WP), PARAMETER :: ONE = 1.0_WP
! Intrinsic functions
      INTRINSIC ABS, LOG

! Set RINF to largest positive real number (infinity)
      RINF = HUGE (RINF)
      RINF = RINF / N
      WARN1 = 0
      WARN2 = 0
      WARN4 = 0
! Check value of JOB
      IF (ICNTL(4).EQ.0) THEN
! Check input data
        IF (JOB.LT.1.OR.JOB.GT.5) THEN
           INFO (1) = - 1
           INFO (2) = JOB
           IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'JOB', JOB
           GOTO 99
        ENDIF
! Check value of N
        IF (N.LT.1) THEN
           INFO (1) = - 2
           INFO (2) = N
           IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'N', N
           GOTO 99
        ENDIF
! Check value of NE
        IF (NE.LT.1) THEN
           INFO (1) = - 3
           INFO (2) = NE
           IF (ICNTL(1) .GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'NE', NE
           GOTO 99
        ENDIF
! Check value of M
        IF (M.LT.N) THEN
           INFO (1) = - 4
           INFO (2) = M
           IF (ICNTL(1) .GE.0 .AND. ICNTL(5).GT.0)   &
             WRITE (ICNTL(1) , 9001) INFO (1) , 'M', M
           GOTO 99
        ENDIF
! Check LIW
!     IF (JOB.EQ.1) K = 4*N +   M
!     IF (JOB.EQ.2) K = 2*N + 2*M
!     IF (JOB.EQ.3) K = 8*N + 2*M + NE
!     IF (JOB.EQ.4) K = 3*N + 2*M
!     IF (JOB.EQ.5) K = 3*N + 2*M
!     IF (JOB.EQ.6) K = 3*N + 2*M + NE
!     IF (LIW.LT.K) THEN
!       INFO(1) = -5
!       INFO(2) = K
!       IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
!       GO TO 99
!     ENDIF
! Check LDW; If JOB = 1, do not check
!     IF (JOB.GT.1) THEN
!       IF (JOB.EQ.2) K =       M
!       IF (JOB.EQ.3) K =           NE
!       IF (JOB.EQ.4) K =     2*M + NE
!       IF (JOB.EQ.5) K = N + 2*M + NE
!       IF (JOB.EQ.6) K = N + 3*M + NE
!       IF (LDW.LT.K) THEN
!         INFO(1) = -5
!         INFO(2) = K
!         IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
!         GO TO 99
!       ENDIF
!     ENDIF
      ENDIF
      IF (ICNTL(4) .EQ.0) THEN
! Check row indices. Use IW(1:M) as workspace
         DO 3 I = 1, M
            IW (I) = 0
    3    END DO
         DO 6 J = 1, N
            DO 4 K = IP (J), IP (J + 1) - 1
               I = IRN (K)
! Check for row indices that are out of range
               IF (I.LT.1.OR.I.GT.M) THEN
                  INFO (1) = - 6
                  INFO (2) = K
                  IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
                    WRITE (ICNTL(1), 9006) INFO (1), J, I
                  GOTO 99
               ENDIF
! Check for repeated row indices within a column
               IF (IW (I) .EQ.J) THEN
                  INFO (1) = - 7
                  INFO (2) = K
                  IF (ICNTL(1).GE.0 .AND. ICNTL(5).GT.0)   &
                    WRITE (ICNTL(1), 9007) INFO (1), J, I
                  GOTO 99
               ELSE
                  IW (I) = J
               ENDIF
    4       END DO
    6    END DO
      ENDIF

! Print diagnostics on input
      IF (ICNTL(5).GT.2) THEN
        IF (ICNTL(3) .GE.0) THEN
          WRITE (ICNTL(3), 9020) JOB, M, N, NE
          IF (ICNTL(5).EQ.3) THEN
            WRITE (ICNTL(3), 9021) (IP (J), J = 1, MIN (10, N + 1) )
            WRITE (ICNTL(3), 9022) (IRN (J), J = 1, MIN (10, NE) )
            IF (JOB.GT.1) WRITE (ICNTL(3), 9023) (A (J), J = 1, MIN &
            (10, NE) )
          ELSE
            WRITE (ICNTL(3), 9021) (IP (J), J = 1, N + 1)
            WRITE (ICNTL(3), 9022) (IRN (J), J = 1, NE)
            IF (JOB.GT.1) WRITE (ICNTL(3), 9023) (A (J), J = 1, NE)
          ENDIF
          WRITE (ICNTL(3), 9024) (ICNTL(J), J = 1, NICNTL)
          WRITE (ICNTL(3), 9025) CNTL(1)
        ENDIF
      ENDIF

! Set components of INFO to zero
      DO I = 1, NINFO
         INFO (I) = 0
      END DO

! Compute maximum matching
      IF (JOB.EQ.1) THEN
! Put length of column J in IW(J)
         DO J = 1, N
            IW (J) = IP (J + 1) - IP (J)
         END DO
! IW(N+1:3N+M+N) is workspace
         CALL MC64_HSL_Z (M, N, IRN, NE, IP, IW (1), PERM, NUM, IW (N +&
         1), IW (2 * N + 1), IW (3 * N + 1), IW (3 * N + M + 1) )
         GOTO 90
      ENDIF

! Compute bottleneck matching
      IF (JOB.EQ.2) THEN
! Pass CNTL(1) to MC64B through DW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
! IW(1:2N+2M), DW(1:M) are workspaces
         CALL MC64_HSL_B (M, N, NE, IP, IRN, A, PERM, NUM, IW (1),     &
         IW (N + 1), IW (2 * N + 1), IW (2 * N + M + 1), DW, RINF)
         GOTO 90
      ENDIF

! Compute bottleneck matching
      IF (JOB.EQ.3) THEN
! Copy IRN(K) into IW(K), ABS(A(K)) into DW(K), K=1..NE
         DO K = 1, NE
            IW (K) = IRN (K)
            DW (K) = ABS (A (K) )
         END DO
! Sort entries in each column by decreasing value.
         CALL MC64RD (N, NE, IP, IW, DW)
! Pass CNTL(1) to MC64S through FACT
         FACT = MAX (ZERO, CNTL (1) )
! IW(NE+1:NE+5N+M+(3N+M)) is workspace
         CALL MC64_HSL_S (M, N, NE, IP, IW (1), DW, PERM, NUM, IW (    &
         NE+1), IW (NE+N + 1), IW (NE+2 * N + 1), IW (NE+3 * N + 1),    &
         IW (NE+4 * N + 1), IW (NE+5 * N + 1), IW (NE+5 * N + M + 1),   &
         FACT, RINF)
         GOTO 90
      ENDIF

      IF (JOB.EQ.4) THEN
         DO J = 1, N
            FACT = ZERO
            DO K = IP (J), IP (J + 1) - 1
               IF (ABS (A (K) ) .GT.FACT) FACT = ABS (A (K) )
            END DO
            DO K = IP (J), IP (J + 1) - 1
               DW (2 * M + K) = FACT - ABS (A (K) )
            END DO
         END DO
! B = DW(2M+1:2M+NE); IW(1:3N+2M) and DW(1:2M) are workspaces
! Pass CNTL(1) to MC64W through DW(1)
! Pass JOB to MC64W through IW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
         IW (1) = JOB
! Call MC64W
         CALL MC64_HSL_W (M, N, NE, IP, IRN, DW (2 * M + 1), PERM, NUM,&
         IW (1), IW (N + 1), IW (2 * N + 1), IW (3 * N + 1), IW (3 * N +&
         M + 1), DW (1), DW (M + 1), RINF)
         GOTO 90
      ENDIF

      IF (JOB.EQ.5.or.JOB.EQ.6) THEN
         IF (JOB.EQ.5) THEN
            DO 75 J = 1, N
               FACT = ZERO
               DO K = IP (J), IP (J + 1) - 1
                  DW (2 * M + N + K) = ABS (A (K) )
                  IF (DW (2 * M + N + K) .GT.FACT) FACT = DW (2 * M + N &
                  + K)
               END DO
               DW (2 * M + J) = FACT
!CC Significant change made here so that column with single
!   zero gets set to RINF and not 1.0
               IF (FACT.NE.ZERO) THEN
                  FACT = LOG (FACT)
               ELSE
                  FACT = RINF
               ENDIF
               DO K = IP (J), IP (J + 1) - 1
                  IF (DW (2 * M + N + K) .NE.ZERO) THEN
                     DW (2 * M + N + K) = FACT - LOG (DW (2 * M + N + K)&
                     )
                  ELSE
!                  write(*,*) 'set diag to ',RINF
                     DW (2 * M + N + K) = RINF
!5.0D+14
!*RINF/(N+1)
                  ENDIF
               END DO
!           ELSE
!             DO 71 K = IP(J),IP(J+1)-1
!               DW(2*M+N+K) = ONE
!  71         CONTINUE
!           ENDIF
   75       END DO
         ENDIF

!       IF (JOB.EQ.6) THEN
!         DO 175 K = 1,NE
!           IW(3*N+2*M+K) = IRN(K)
!           DW(2*M+N+K) = ABS(A(K))
! 175     CONTINUE
!         DO 61 I = 1,M
!           DW(2*M+N+NE+I) = ZERO
!  61     CONTINUE
!         DO 63 J = 1,N
!           DO 62 K = IP(J),IP(J+1)-1
!             I = IRN(K)
!             IF (DW(2*M+N+K).GT.DW(2*M+N+NE+I)) THEN
!               DW(2*M+N+NE+I) = DW(2*M+N+K)
!             ENDIF
!  62       CONTINUE
!  63     CONTINUE
!         DO 64 I = 1,M
!           IF (DW(2*M+N+NE+I).NE.ZERO) THEN
!             DW(2*M+N+NE+I) = 1/DW(2*M+N+NE+I)
!           ENDIF
!  64     CONTINUE
!         DO 66 J = 1,N
!           DO 65 K = IP(J),IP(J+1)-1
!             I = IRN(K)
!             DW(2*M+N+K) = DW(2*M+N+NE+I) * DW(2*M+N+K)
!  65       CONTINUE
!  66     CONTINUE
!         CALL MC64R(N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1))
!         DO 176 J = 1,N
!           IF (IP(J).NE.IP(J+1)) THEN
!             FACT = DW(2*M+N+IP(J))
!           ELSE
!             FACT = ZERO
!           ENDIF
!           DW(2*M+J) = FACT
!           IF (FACT.NE.ZERO) THEN
!             FACT = LOG(FACT)
!             DO 170 K = IP(J),IP(J+1)-1
!               IF (DW(2*M+N+K).NE.ZERO) THEN
!                 DW(2*M+N+K) = FACT - LOG(DW(2*M+N+K))
!               ELSE
!                  write(*,*) 'set diag to ',RINF
!                 DW(2*M+N+K) = RINF
!5.0D+14
!0.5* RINF/(N+1)
!               ENDIF
! 170         CONTINUE
!           ELSE
!             DO 171 K = IP(J),IP(J+1)-1
!               DW(2*M+N+K) = ONE
! 171         CONTINUE
!           ENDIF
! 176     CONTINUE
!       ENDIF
! Pass CNTL(1) to MC64W through DW(1)
! Pass JOB to MC64W through IW(1)
         DW (1) = MAX (ZERO, CNTL (1) )
         IW (1) = JOB
! Call MC64W
         IF (JOB.EQ.5) THEN
            CALL MC64_HSL_W (M, N, NE, IP, IRN, DW (2 * M + N + 1),    &
            PERM, NUM, IW (1), IW (N + 1), IW (2 * N + 1), IW (3 * N +  &
            1), IW (3 * N + M + 1), DW (1), DW (M + 1), RINF)
         ENDIF

!       IF (JOB.EQ.6) THEN
!         CALL MC64_HSL_W(M,N,NE,IP,IW(3*N+2*M+1),DW(2*M+N+1),PERM,NUM,
!    &         IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(3*N+M+1),
!    &         DW(1),DW(M+1),RINF)
!       ENDIF
!       IF (JOB.EQ.6) THEN
!         DO 79 I = 1,M
!           IF (DW(2*M+N+NE+I).NE.0) THEN
!             DW(I) = DW(I) + LOG(DW(2*M+N+NE+I))
!           ENDIF
!  79     CONTINUE
!       ENDIF
         IF (NUM.EQ.N) THEN
            DO    J = 1, N
               IF (DW (2 * M + J) .NE.ZERO) THEN
                  DW (M + J) = DW (M + J) - LOG (DW (2 * M + J) )
               ELSE
                  DW (M + J) = ZERO
               ENDIF
            END DO
         ENDIF
! Check size of row and column scaling factors
         FACT = 0.5 * LOG (RINF)
         DO J = 1, N
            IF (DW (M + J) .LT.FACT) CYCLE
            WARN2 = 2
! Scaling factor is large, return with warning
            GOTO 90
         END DO
         DO I = 1, M
            IF (DW (I) .LT.FACT) CYCLE
            WARN2 = 2
! Scaling factor is large, return with warning
            GOTO 90
         END DO
!       GO TO 90
      ENDIF

! If matrix is structurally singular, return with warning
   90 IF (NUM.LT.N) WARN1 = 1

! If CNTL(1) is negative and treated as zero, return with warning
      IF (JOB.EQ.4.OR.JOB.EQ.5.OR.JOB.EQ.6) THEN
         IF (CNTL (1) .LT.ZERO) WARN4 = 4
      ENDIF

! Set warning flag and print warnings (only if no errors were found)
      IF (INFO (1) .EQ.0) THEN
         INFO (1) = WARN1 + WARN2 + WARN4
         IF (INFO (1).GT.0 .AND. ICNTL(2).GE.0 .AND. ICNTL(5).GT.1) THEN
           WRITE (ICNTL(2), 9010) INFO (1)
           IF (WARN1.EQ.1) WRITE (ICNTL(2), 9011)
           IF (WARN2.EQ.2) WRITE (ICNTL(2), 9012)
           IF (WARN4.EQ.4) WRITE (ICNTL(2), 9014)
         ENDIF
      ENDIF
! Print diagnostics on output
      IF (ICNTL(5).GT.2) THEN
        IF (ICNTL(3).GE.0) THEN
          WRITE (ICNTL(3), 9030) (INFO (J), J = 1, 2)
          WRITE (ICNTL(3), 9031) NUM
          IF (ICNTL(5) .EQ.3) THEN
            WRITE (ICNTL(3), 9032) (PERM (J), J = 1, MIN (10, M) )
            IF (JOB.EQ.5.OR.JOB.EQ.6) THEN
              WRITE (ICNTL(3), 9033) (DW (J), J = 1, MIN (10, M) )
              WRITE (ICNTL(3), 9034) (DW (M + J), J = 1, MIN (10, N))
            ENDIF
          ELSE
            WRITE (ICNTL(3), 9032) (PERM (J), J = 1, M)
            IF (JOB.EQ.5.OR.JOB.EQ.6) THEN
              WRITE (ICNTL(3), 9033) (DW (J), J = 1, M)
              WRITE (ICNTL(3), 9034) (DW (M + J), J = 1, N)
            ENDIF
          ENDIF
        ENDIF
      ENDIF

! Return from subroutine.
   99 RETURN

 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,               &
     &        ' because ',(A),' = ',I10)
!9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
!    &        '        LIW too small, must be at least ',I8)
!9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
!    &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
     &        '        Column ',I8,                                     &
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/               &
     &        '        Column ',I8,                                     &
     &        ' contains two or more entries with row index ',I8)

 9010 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2)
 9011 FORMAT ('        - The matrix is structurally singular.')
 9012 FORMAT ('        - Some scaling factors may be too large.')
 9014 FORMAT ('        - CNTL(1) is negative and was treated as zero.')

 9020 FORMAT (' ****** Input parameters for MC64AD:'/                   &
     &        ' JOB =',I10/' M   =',I10/' N   =',I10/' NE  =',I10)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9024 FORMAT (' ICNTL(1:10)= ',8I8/(14X,2I8))
 9025 FORMAT (' CNTL(1)    = ',1PD14.4)
 9030 FORMAT (' ****** Output parameters for MC64AD:'/                  &
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' PERM(1:M)  = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:M)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(M+1:M+N)= ',5(F11.3)/(14X,5(F11.3)))
END SUBROUTINE MC64_HSL_A

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_B (M, N, NE, IP, IRN, A, IPERM, NUM, JPERM,  &
      PR, Q, L, D, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUM
      INTEGER IP (N + 1), IRN (NE), IPERM (M), JPERM (N), PR (N),       &
      Q (M), L (M)
      REAL(WP) A (NE)
      REAL(WP) D (M), RINF

! N, NE, IP, IRN are described in MC64AD.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry
!   that corresponds to IRN(K). It is not altered.
! IPERM is an INTEGER array of length M. On exit, it contains the
!    matching: IPERM(I) = 0 or row I is matched to column IPERM(I).
! NUM is INTEGER variable. On exit, it contains the cardinality of the
!    matching stored in IPERM.
! D is a REAL work array of length M.
!    On entry, D(1) contains the relaxation parameter RLX.
! RINF is the largest positive real number

! Local variables
      INTEGER I, II, J, JJ, JORD, Q0, QLEN, IDUM, JDUM, ISP, JSP, K, KK,&
      KK1, KK2, I0, UP, LOW, LPOS
      REAL(WP) CSP, DI, DNEW, DQ0, AI, A0, BV, TBV, RLX
! Local parameters
      REAL(WP), PARAMETER :: ZERO=0.0_WP
      REAL(WP), PARAMETER :: MINONE=-1.0_WP
      INTRINSIC ABS, MIN
! External subroutines and/or functions
      EXTERNAL MC64DD, MC64ED, MC64FD

! Initialize variables to eliminate copmiler warnings
      I0 = -HUGE(I0); ISP = -HUGE(ISP); JSP = -HUGE(JSP)

      RLX = D (1)
! Initialization
      NUM = 0
      BV = RINF
      DO K = 1, N
         JPERM (K) = 0
         PR (K) = IP (K)
      END DO
      DO K = 1, M
         IPERM (K) = 0
         D (K) = ZERO
      END DO

      DO J = 1, N
         A0 = MINONE
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            AI = ABS (A (K) )
            IF (AI.GT.D (I) ) D (I) = AI
            IF (JPERM (J) .NE.0) CYCLE
            IF (AI.GE.BV) THEN
               A0 = BV
               IF (IPERM (I) .NE.0) CYCLE
               JPERM (J) = I
               IPERM (I) = J
               NUM = NUM + 1
            ELSE
               IF (AI.LE.A0) CYCLE
               A0 = AI
               I0 = I
            ENDIF
         END DO
         IF (A0.NE.MINONE.AND.A0.LT.BV) THEN
            BV = A0
            IF (IPERM (I0) .NE.0) CYCLE
            IPERM (I0) = J
            JPERM (J) = I0
            NUM = NUM + 1
         ENDIF
      END DO

      IF (M.EQ.N) THEN
! Update BV with smallest of all the largest maximum absolute values
! of the rows. D(I) contains the largest absolute value in row I.
         DO I = 1, M
            BV = MIN (BV, D (I) )
         END DO
      ENDIF

! Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GOTO 1000

! Rescan unassigned columns; improve initial assignment
      DO J = 1, N
         IF (JPERM (J) .NE.0) CYCLE
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            AI = ABS (A (K) )
            IF (AI.LT.BV) CYCLE
            IF (IPERM (I) .EQ.0) GOTO 90
            JJ = IPERM (I)
            KK1 = PR (JJ)
            KK2 = IP (JJ + 1) - 1
            IF (KK1.GT.KK2) CYCLE
            DO KK = KK1, KK2
               II = IRN (KK)
               IF (IPERM (II) .NE.0) CYCLE
               IF (ABS (A (KK) ) .GE.BV) GOTO 80
            END DO
            PR (JJ) = KK2 + 1
         END DO
         CYCLE
   80    JPERM (JJ) = II
         IPERM (II) = JJ
         PR (JJ) = KK + 1
   90    NUM = NUM + 1
         JPERM (J) = I
         IPERM (I) = J
         PR (J) = K + 1
      END DO

! Shortcut if all columns are matched at this stage.
      IF (NUM.EQ.N) GOTO 1000

! Prepare for main loop
      DO I = 1, M
         D (I) = MINONE
         L (I) = 0
      END DO
! TBV is a relaxed value of BV (ie TBV is slightly smaller than BV).
      TBV = BV * (1 - RLX)

! Main loop ... each pass round this loop is similar to Dijkstra's
! algorithm for solving the single source shortest path problem

      DO JORD = 1, N

         IF (JPERM (JORD) .NE.0) CYCLE
         QLEN = 0
         LOW = M + 1
         UP = M + 1
! CSP is cost of shortest path to any unassigned row
! ISP is matrix position of unassigned row element in shortest path
! JSP is column index of unassigned row element in shortest path
         CSP = MINONE
! Build shortest path tree starting from unassigned column JORD
         J = JORD
         PR (J) = - 1

! Scan column J
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            DNEW = ABS (A (K) )
            IF (CSP.GE.DNEW) CYCLE
            IF (IPERM (I) .EQ.0) THEN
! Row I is unassigned; update shortest path info
               CSP = DNEW
               ISP = I
               JSP = J
               IF (CSP.GE.TBV) GOTO 160
            ELSE
               D (I) = DNEW
               IF (DNEW.GE.TBV) THEN
! Add row I to Q2
                  LOW = LOW - 1
                  Q (LOW) = I
               ELSE
! Add row I to Q, and push it
                  QLEN = QLEN + 1
                  L (I) = QLEN
                  CALL MC64DD (I, M, Q, D, L, 1)
               ENDIF
               JJ = IPERM (I)
               PR (JJ) = J
            ENDIF
         END DO

         DO JDUM = 1, NUM
! If Q2 is empty, extract new rows from Q
            IF (LOW.EQ.UP) THEN
               IF (QLEN.EQ.0) EXIT
               I = Q (1)
               IF (CSP.GE.D (I) ) EXIT
               BV = D (I)
               TBV = BV * (1 - RLX)
               DO IDUM = 1, M
                  CALL MC64ED (QLEN, M, Q, D, L, 1)
                  L (I) = 0
                  LOW = LOW - 1
                  Q (LOW) = I
                  IF (QLEN.EQ.0) EXIT
                  I = Q (1)
                  IF (D (I) .LT.TBV) EXIT
               END DO
! End of dummy loop; this point is never reached
            ENDIF
! Move row Q0
            UP = UP - 1
            Q0 = Q (UP)
            DQ0 = D (Q0)
            L (Q0) = UP
! Scan column that matches with row Q0
            J = IPERM (Q0)
            DO K = IP (J), IP (J + 1) - 1
               I = IRN (K)
! Update D(I); only if row I is not marked
               IF (L (I) .GE.UP) CYCLE
               DNEW = MIN (DQ0, ABS (A (K) ) )
               IF (CSP.GE.DNEW) CYCLE
               IF (IPERM (I) .EQ.0) THEN
! Row I is unassigned; update shortest path info
                  CSP = DNEW
                  ISP = I
                  JSP = J
                  IF (CSP.GE.TBV) GOTO 160
               ELSE
                  DI = D (I)
                  IF (DI.GE.TBV.OR.DI.GE.DNEW) CYCLE
                  D (I) = DNEW
                  IF (DNEW.GE.TBV) THEN
! Delete row I from Q (if necessary); add row I to Q2
                     IF (DI.NE.MINONE) THEN
                        LPOS = L (I)
                        CALL MC64FD (LPOS, QLEN, M, Q, D, L, 1)
                     ENDIF
                     L (I) = 0
                     LOW = LOW - 1
                     Q (LOW) = I
                  ELSE
! Add row I to Q (if necessary); push row I up Q
                     IF (DI.EQ.MINONE) THEN
                        QLEN = QLEN + 1
                        L (I) = QLEN
                     ENDIF
                     CALL MC64DD (I, M, Q, D, L, 1)
                  ENDIF
! Update tree
                  JJ = IPERM (I)
                  PR (JJ) = J
               ENDIF
            END DO
         END DO

! If CSP = MINONE, no augmenting path is found
  160    IF (CSP.EQ.MINONE) GOTO 190
! Update bottleneck value
         BV = MIN (BV, CSP)
         TBV = BV * (1 - RLX)
! Find augmenting path by tracing backward in PR; update IPERM,JPERM
         NUM = NUM + 1
         I = ISP
         J = JSP
         DO JDUM = 1, NUM + 1
            I0 = JPERM (J)
            JPERM (J) = I
            IPERM (I) = J
            J = PR (J)
            IF (J.EQ. - 1) EXIT
            I = I0
         END DO
! End of dummy loop; this point is never reached
  190    DO 191 KK = UP, M
            I = Q (KK)
            D (I) = MINONE
            L (I) = 0
  191    END DO
         DO 192 KK = LOW, UP - 1
            I = Q (KK)
            D (I) = MINONE
  192    END DO
         DO 193 KK = 1, QLEN
            I = Q (KK)
            D (I) = MINONE
            L (I) = 0
  193    END DO

      END DO
! End of main loop

! BV is now bottleneck value of final matching

! IPERM is complete if M = N and NUM = N
 1000 IF (M.EQ.N.and.NUM.EQ.N) GOTO 2000

! Complete IPERM; L, JPERM are work arrays
      CALL MC64_HSL_X (M, N, IPERM, L, JPERM)

 2000 RETURN
END SUBROUTINE MC64_HSL_B

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_S (M, N, NE, IP, IRN, A, IPERM, NUMX, W, LEN,&
      LENL, LENH, FC, IW, IW4, RLX, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUMX
      INTEGER IP (N + 1), IRN (NE), IPERM (M), W (N), LEN (N), LENL (N),&
      LENH (N), FC (N), IW (M), IW4 (3 * N + M)
      REAL(WP) A (NE), RLX, RINF

! M, N, NE, IP, IRN, are described in MC64AD.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry that
!   corresponds to IRN(k). The entries in each column must be
!   non-negative and ordered by decreasing value.
! IPERM is an INTEGER array of length M. On exit, it contains the
!   bottleneck matching: IPERM(I) - 0 or row I is matched to column
!   IPERM(I).
! NUMX is an INTEGER variable. On exit, it contains the cardinality
!   of the matching stored in IPERM.

! FC is an integer array of length N that contains the list of
!   unmatched columns.
! LEN(J), LENL(J), LENH(J) are integer arrays of length N that point
!   to entries in matrix column J.
!   In the matrix defined by the column parts IP(J)+LENL(J) we know
!   a matching does not exist; in the matrix defined by the column
!   parts IP(J)+LENH(J) we know one exists.
!   LEN(J) lies between LENL(J) and LENH(J) and determines the matrix
!   that is tested for a maximum matching.
! W is an integer array of length N and contains the indices of the
!   columns for which LENL /= LENH.
! WLEN is number of indices stored in array W.
! IW is integer work array of length M.
! IW4 is integer work array of length 3N+M used by MC64U.
!
! RLX is a REAL variable. It is a relaxation
!   parameter for finding the optimal matching.
!
! RINF is the largest positive real number

      INTEGER NUM, NVAL, WLEN, II, I, J, K, L, CNT, MOD, IDUM1, IDUM2,  &
      IDUM3
      REAL(WP) BVAL, BMIN, BMAX
! External subroutines and/or functions
      EXTERNAL MC64QD
! Intrinsic functions

! BMIN and BMAX are such that a maximum matching exists for the input
!   matrix in which all entries smaller than BMIN are dropped.
!   For BMAX, a maximum matching does not exist.
! BVAL is a value between BMIN and BMAX.
! CNT is the number of calls made to MC64U so far.
! NUM is the cardinality of last matching found.

! Compute a first maximum matching from scratch on whole matrix.
      DO J = 1, N
         FC (J) = J
         LEN (J) = IP (J + 1) - IP (J)
      END DO
      DO I = 1, M
         IW (I) = 0
      END DO
! The first call to MC64U
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64_HSL_U (CNT, MOD, M, N, IRN, NE, IP, LEN, FC, IW, NUMX, &
      N, IW4 (1), IW4 (N + 1), IW4 (2 * N + 1), IW4 (2 * N + M + 1) )

! IW contains a maximum matching of length NUMX.
      NUM = NUMX

      IF (NUM.NE.N) THEN
! Matrix is structurally singular
         BMAX = RINF
      ELSE
! Matrix is structurally nonsingular, NUM=NUMX=N;
! Set BMAX just above the smallest of all the maximum absolute
! values of the columns
         BMAX = RINF
         DO J = 1, N
            BVAL = 0.0
            DO K = IP (J), IP (J + 1) - 1
               IF (A (K) .GT.BVAL) BVAL = A (K)
            END DO
            IF (BVAL.LT.BMAX) BMAX = BVAL
         END DO
! ... should print warning if BMAX == RINF
         BMAX = 1.001 * BMAX
      ENDIF

! Initialize BVAL,BMIN
      BVAL = 0.0
      BMIN = 0.0
! Initialize LENL,LEN,LENH,W,WLEN according to BMAX.
! Set LEN(J), LENH(J) just after last entry in column J.
! Set LENL(J) just after last entry in column J with value >= BMAX.
      WLEN = 0
      DO J = 1, N
         L = IP (J + 1) - IP (J)
         LENH (J) = L
         LEN (J) = L
         DO K = IP (J), IP (J + 1) - 1
            IF (A (K) .LT.BMAX) GOTO 46
         END DO
! Column J is empty or all entries are >= BMAX
         K = IP (J + 1)
   46    LENL (J) = K - IP (J)
! Add J to W if LENL(J) /= LENH(J)
         IF (LENL (J) .EQ.L) CYCLE
         WLEN = WLEN + 1
         W (WLEN) = J
      END DO

! Main loop
      DO IDUM1 = 1, NE
         IF (NUM.EQ.NUMX) THEN
! We have a maximum matching in IW; store IW in IPERM
            DO  I = 1, M
               IPERM (I) = IW (I)
            END DO
! Keep going round this loop until matching IW is no longer maximum.
            DO 80 IDUM2 = 1, NE
               BMIN = BVAL
               IF (BMAX - BMIN.LE.RLX) GOTO 1000
! Find splitting value BVAL
               CALL MC64QD (IP, LENL, LEN, W, WLEN, A, NVAL, BVAL)
               IF (NVAL.LE.1) GOTO 1000
! Set LEN such that all matrix entries with value < BVAL are
! discarded. Store old LEN in LENH. Do this for all columns W(K).
! Each step, either K is incremented or WLEN is decremented.
               K = 1
               DO IDUM3 = 1, N
                  IF (K.GT.WLEN) EXIT
                  J = W (K)
                  DO II = IP (J) + LEN (J) - 1, IP (J) + LENL (J),   &
                  - 1
                     IF (A (II) .GE.BVAL) EXIT
                     I = IRN (II)
                     IF (IW (I) .NE.J) CYCLE
! Remove entry from matching
                     IW (I) = 0
                     NUM = NUM - 1
                     FC (N - NUM) = J
                  END DO
                  LENH (J) = LEN (J)
! IP(J)+LEN(J)-1 is last entry in column >= BVAL
                  LEN (J) = II - IP (J) + 1
! If LENH(J) = LENL(J), remove J from W
                  IF (LENL (J) .EQ.LENH (J) ) THEN
                     W (K) = W (WLEN)
                     WLEN = WLEN - 1
                  ELSE
                     K = K + 1
                  ENDIF
               END DO
               IF (NUM.LT.NUMX) EXIT
   80       END DO
! End of dummy loop; this point is never reached
! Set mode for next call to MC64U
            MOD = 1
         ELSE
! We do not have a maximum matching in IW.
            BMAX = BVAL
! BMIN is the bottleneck value of a maximum matching;
! for BMAX the matching is not maximum, so BMAX>BMIN
! and following condition is always false if RLX = 0.0
            IF (BMAX - BMIN.LE.RLX) GOTO 1000
! Find splitting value BVAL
            CALL MC64QD (IP, LEN, LENH, W, WLEN, A, NVAL, BVAL)
            IF (NVAL.EQ.0.OR.BVAL.EQ.BMIN) GOTO 1000
! Set LEN such that all matrix entries with value >= BVAL are
! inside matrix. Store old LEN in LENL. Do this for all columns W(K).
! Each step, either K is incremented or WLEN is decremented.
            K = 1
            DO 87 IDUM3 = 1, N
               IF (K.GT.WLEN) GOTO 88
               J = W (K)
               DO 85 II = IP (J) + LEN (J), IP (J) + LENH (J) - 1
                  IF (A (II) .LT.BVAL) GOTO 86
   85          END DO
   86          LENL (J) = LEN (J)
               LEN (J) = II - IP (J)
               IF (LENL (J) .EQ.LENH (J) ) THEN
                  W (K) = W (WLEN)
                  WLEN = WLEN - 1
               ELSE
                  K = K + 1
               ENDIF
   87       END DO
! End of dummy loop; this point is never reached
! Set mode for next call to MC64U
   88       MOD = 0
         ENDIF
         CNT = CNT + 1
         CALL MC64_HSL_U (CNT, MOD, M, N, IRN, NE, IP, LEN, FC, IW,    &
         NUM, NUMX, IW4 (1), IW4 (N + 1), IW4 (2 * N + 1), IW4 (2 * N + &
         M + 1) )

! IW contains maximum matching of length NUM
      END DO
! End of dummy loop; this point is never reached

! BMIN is now bottleneck value of final matching

! IPERM is complete if M = N and NUMX = N
 1000 IF (M.EQ.N.and.NUMX.EQ.N) GOTO 2000

! Complete IPERM; IW, W are work arrays
      CALL MC64_HSL_X (M, N, IPERM, IW, W)

 2000 RETURN
END SUBROUTINE MC64_HSL_S

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_U (ID, MOD, M, N, IRN, LIRN, IP, LENC, FC,   &
      IPERM, NUM, NUMX, PR, ARP, CV, OUT)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER ID, MOD, M, N, LIRN, NUM, NUMX
      INTEGER ARP (N), CV (M), IRN (LIRN), IP (N), FC (N), IPERM (M),   &
      LENC (N), OUT (N), PR (N)

! PR(J) is the previous column to J in the depth first search.
!   Array PR is used as workspace in the sorting algorithm.
! Elements (I,IPERM(I)) I=1,..,M are entries at the end of the
!   algorithm unless N assignments have not been made in which case
!   N-NUM pairs (I,IPERM(I)) will not be entries in the matrix.
! CV(I) is the most recent loop number (ID+JORD) at which row I
!   was visited.
! ARP(J) is the number of entries in column J which have been scanned
!   when looking for a cheap assignment.
! OUT(J) is one less than the number of entries in column J which have
!   not been scanned during one pass through the main loop.
! NUMX is maximum possible size of matching.

      INTEGER I, II, IN1, IN2, J, J1, JORD, K, KK, LAST, NFC, NUM0,     &
      NUM1, NUM2, ID0, ID1

! Initialize variables to eliminate copmiler warnings
      I = -1; II = -1

      IF (ID.EQ.1) THEN
! The first call to MC64U.
! Initialize CV and ARP; parameters MOD, NUMX are not accessed
         DO I = 1, M
            CV (I) = 0
         END DO
         DO J = 1, N
            ARP (J) = 0
         END DO
         NUM1 = N
         NUM2 = N
      ELSE
! Not the first call to MC64U.
! Re-initialize ARP if entries were deleted since last call to MC64U
         IF (MOD.EQ.1) THEN
            DO J = 1, N
               ARP (J) = 0
            END DO
         ENDIF
         NUM1 = NUMX
         NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM

! NUM0 is size of input matching
! NUM1 is maximum possible size of matching
! NUM2 is maximum allowed number of unassigned rows/columns
! NUM is size of current matching

! Quick return if possible
!      IF (NUM.EQ.N) GO TO 199
! NFC is number of rows/columns that could not be assigned
      NFC = 0
! Integers ID0+1 to ID0+N are unique numbers for call ID to MC64U,
! so 1st call uses 1..N, 2nd call uses N+1..2N, etc
      ID0 = (ID-1) * N

! Main loop. Each pass round this loop either results in a new
! assignment or gives a column with no assignment

      OUTER: DO JORD = NUM0 + 1, N

! Each pass uses unique number ID1
         ID1 = ID0 + JORD
! J is unmatched column
         J = FC (JORD-NUM0)
         PR (J) = - 1
         SCAN: DO K = 1, JORD
! Look for a cheap assignment
            IF (ARP (J) .GE.LENC (J) ) GOTO 30
            IN1 = IP (J) + ARP (J)
            IN2 = IP (J) + LENC (J) - 1
            DO II = IN1, IN2
               I = IRN (II)
               IF (IPERM (I) .EQ.0) GOTO 80
            END DO
! No cheap assignment in row
            ARP (J) = LENC (J)
! Begin looking for assignment chain starting with row J
   30       OUT (J) = LENC (J) - 1
! Inner loop.  Extends chain by one or backtracks
            DO KK = 1, JORD
               IN1 = OUT (J)
               IF (IN1.LT.0) GOTO 50
               IN2 = IP (J) + LENC (J) - 1
               IN1 = IN2 - IN1
! Forward scan
               DO II = IN1, IN2
                  I = IRN (II)
                  IF (CV (I) .EQ.ID1) CYCLE
! Column J has not yet been accessed during this pass
                  J1 = J
                  J = IPERM (I)
                  CV (I) = ID1
                  PR (J) = J1
                  OUT (J1) = IN2 - II - 1
                  CYCLE SCAN
               END DO
! Backtracking step.
   50          J1 = PR (J)
               IF (J1.EQ. - 1) THEN
! No augmenting path exists for column J.
                  NFC = NFC + 1
                  FC (NFC) = J
                  IF (NFC.GT.NUM2) THEN
! A matching of maximum size NUM1 is not possible
                     LAST = JORD
                     GOTO 101
                  ENDIF
                  CYCLE OUTER
               ENDIF
               J = J1
            END DO
! End of dummy loop; this point is never reached
         END DO SCAN
! End of dummy loop; this point is never reached

! New assignment is made.
   80    IPERM (I) = J
         ARP (J) = II - IP (J) + 1
         NUM = NUM + 1
         DO 90 K = 1, JORD
            J = PR (J)
            IF (J.EQ. - 1) GOTO 95
            II = IP (J) + LENC (J) - OUT (J) - 2
            I = IRN (II)
            IPERM (I) = J
   90    END DO
! End of dummy loop; this point is never reached

   95    IF (NUM.EQ.NUM1) THEN
! A matching of maximum size NUM1 is found
            LAST = JORD
            GOTO 101
         ENDIF
!
      END DO OUTER

! All unassigned columns have been considered
      LAST = N

! Now, a transversal is computed or is not possible.
! Complete FC before returning.
  101 DO 110 JORD = LAST + 1, N
         NFC = NFC + 1
         FC (NFC) = FC (JORD-NUM0)
  110 END DO

!  199 RETURN
      RETURN
END SUBROUTINE MC64_HSL_U

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_W (M, N, NE, IP, IRN, A, IPERM, NUM, JPERM,  &
      OUT, PR, Q, L, U, D, RINF)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
      INTEGER M, N, NE, NUM
      INTEGER IP (N + 1), IRN (NE), IPERM (M), JPERM (N), OUT (N),      &
      PR (N), Q (M), L (M)
      REAL(WP) A (NE), U (M), D (M), RINF

! M, N, NE, IP, IRN are described in MC64AD.
! A is a REAL array of length NE.
!   A(K), K=1..NE, must be set to the value of the entry that
!   corresponds to IRN(K). It is not altered.
!   All values A(K) must be non-negative.
! IPERM is an INTEGER array of length M. On exit, it contains the
!   weighted matching: IPERM(I) = 0 or row I is matched to column
!   IPERM(I).
! NUM is an INTEGER variable. On exit, it contains the cardinality of
!   the matching stored in IPERM.
! D is a REAL array of length M.
!   On exit, V = D(1:N) contains the dual column variable.
!   If U(1:M) denotes the dual row variable and if the matrix
!   is structurally nonsingular (NUM = N), the following holds:
!      U(I)+V(J) <= A(I,J)  if IPERM(I) /= J
!      U(I)+V(J)  = A(I,J)  if IPERM(I)  = J
!      U(I) = 0  if IPERM(I) = 0
!      V(J) = 0  if there is no I for which IPERM(I) = J
!   On entry, U(1) contains the relaxation parameter RLX.
! RINF is the largest positive real number

! Local variables
      INTEGER I, I0, II, J, JJ, JORD, Q0, QLEN, JDUM, ISP, JSP, K, K0,  &
      K1, K2, KK, KK1, KK2, UP, LOW, LPOS
      REAL(WP) CSP, DI, DMIN, DNEW, DQ0, VJ
! Local parameters
      REAL(WP),PARAMETER :: ZERO=0.0_WP
! External subroutines and/or functions
      EXTERNAL MC64DD, MC64ED, MC64FD

! Initialize variables to eliminate copmiler warnings
      ISP = -HUGE(ISP); JSP = -HUGE(JSP)

! Set RINF to largest positive real number
      RINF = HUGE(RINF)

! Initialization
      NUM = 0
      DO K = 1, N
         D (K) = ZERO
         JPERM (K) = 0
         PR (K) = IP (K)
      END DO
      DO K = 1, M
         U (K) = RINF
         IPERM (K) = 0
         L (K) = 0
      END DO
! Initialize U(I)
      DO J = 1, N
         DO K = IP (J), IP (J + 1) - 1
            I = IRN (K)
            IF (A (K) .GT.U (I) ) CYCLE
            U (I) = A (K)
            IPERM (I) = J
            L (I) = K
         END DO
      END DO
      DO I = 1, M
         J = IPERM (I)
         IF (J.EQ.0) CYCLE
! Row I is not empty
         IPERM (I) = 0
         IF (JPERM (J) .NE.0) CYCLE
! Don't choose cheap assignment from dense columns
         IF (IP (J + 1) - IP (J) .GT.N / 10.AND.N.GT.50) CYCLE
! Assignment of column J to row I
         NUM = NUM + 1
         IPERM (I) = J
         JPERM (J) = L (I)
      END DO

!      write(14,*) 'Number of cheap assignments ',NUM

      IF (NUM.EQ.N) GOTO 1000
! Scan unassigned columns; improve assignment
      DO J = 1, N
! JPERM(J) ne 0 iff column J is already assigned
         IF (JPERM (J) .NE.0) CYCLE
         K1 = IP (J)
         K2 = IP (J + 1) - 1
! Continue only if column J is not empty
         IF (K1.GT.K2) CYCLE
!       VJ = RINF
! Changes made to allow for NaNs
         I0 = IRN (K1)
         VJ = A (K1) - U (I0)
         K0 = K1
         DO K = K1 + 1, K2
            I = IRN (K)
            DI = A (K) - U (I)
            IF (DI.GT.VJ) CYCLE
            IF (.NOT.(DI.LT.VJ.OR.DI.EQ.RINF)) THEN
               IF (IPERM (I) .NE.0.OR.IPERM (I0) .EQ.0) CYCLE
            ENDIF
            VJ = DI
            I0 = I
            K0 = K
         END DO
         D (J) = VJ
         K = K0
         I = I0
         IF (IPERM (I) .EQ.0) GOTO 90
         DO K = K0, K2
            I = IRN (K)
            IF (A (K) - U (I) .GT.VJ) CYCLE
            JJ = IPERM (I)
! Scan remaining part of assigned column JJ
            KK1 = PR (JJ)
            KK2 = IP (JJ + 1) - 1
            IF (KK1.GT.KK2) CYCLE
            DO KK = KK1, KK2
               II = IRN (KK)
               IF (IPERM (II) .GT.0) CYCLE
               IF (A (KK) - U (II) .LE.D (JJ) ) GOTO 80
            END DO
            PR (JJ) = KK2 + 1
         END DO
         CYCLE
   80    JPERM (JJ) = KK
         IPERM (II) = JJ
         PR (JJ) = KK + 1
   90    NUM = NUM + 1
         JPERM (J) = K
         IPERM (I) = J
         PR (J) = K + 1
      END DO

!     write(14,*) 'Number of improved  assignments ',NUM

      IF (NUM.EQ.N) GOTO 1000

! Prepare for main loop
      DO I = 1, M
         D (I) = RINF
         L (I) = 0
      END DO

! Main loop ... each pass round this loop is similar to Dijkstra's
! algorithm for solving the single source shortest path problem

      DO JORD = 1, N

         IF (JPERM (JORD) .NE.0) CYCLE
! JORD is next unmatched column
! DMIN is the length of shortest path in the tree
         DMIN = RINF
         QLEN = 0
         LOW = N + 1
         UP = N + 1
! CSP is the cost of the shortest augmenting path to unassigned row
! IRN(ISP). The corresponding column index is JSP.
         CSP = RINF
! Build shortest path tree starting from unassigned column (root) JORD
         J = JORD
         PR (J) = - 1

! Scan column J
         DO K = IP (J), IP (J + 1) - 1
!         IF (N.EQ.3) THEN
!           write(14,*) 'Scanning column ',J
!           write(14,*) 'IP  ',IP(1:4)
!           write(14,*) 'IRN ',IRN(1:6)
!           write(14,*) 'A   ',A(1:6)
!           write(14,*) 'U ',U(1:3)
!           write(14,*) 'IPERM ',IPERM(1:3)
!         ENDIF
            I = IRN (K)
            DNEW = A (K) - U (I)
            IF (DNEW.GE.CSP) CYCLE
            IF (IPERM (I) .EQ.0) THEN
               CSP = DNEW
               ISP = K
               JSP = J
            ELSE
               IF (DNEW.LT.DMIN) DMIN = DNEW
               D (I) = DNEW
               QLEN = QLEN + 1
               Q (QLEN) = K
            ENDIF
         END DO
! Initialize heap Q and Q2 with rows held in Q(1:QLEN)
         Q0 = QLEN
         QLEN = 0
         DO KK = 1, Q0
            K = Q (KK)
            I = IRN (K)
            IF (CSP.LE.D (I) ) THEN
               D (I) = RINF
               CYCLE
            ENDIF
            IF (D (I) .LE.DMIN) THEN
               LOW = LOW - 1
               Q (LOW) = I
               L (I) = LOW
            ELSE
               QLEN = QLEN + 1
               L (I) = QLEN
               CALL MC64DD (I, M, Q, D, L, 2)
            ENDIF
! Update tree
            JJ = IPERM (I)
            OUT (JJ) = K
            PR (JJ) = J
         END DO

         DO JDUM = 1, NUM

! If Q2 is empty, extract rows from Q
            IF (LOW.EQ.UP) THEN
               IF (QLEN.EQ.0) GOTO 160
               I = Q (1)
               IF (D (I) .GE.CSP) GOTO 160
               DMIN = D (I)
  152          CALL MC64ED (QLEN, M, Q, D, L, 2)
               LOW = LOW - 1
               Q (LOW) = I
               L (I) = LOW
               IF (QLEN.EQ.0) GOTO 153
               I = Q (1)
               IF (D (I) .GT.DMIN) GOTO 153
               GOTO 152
            ENDIF
! Q0 is row whose distance D(Q0) to the root is smallest
  153       Q0 = Q (UP - 1)
            DQ0 = D (Q0)
! Exit loop if path to Q0 is longer than the shortest augmenting path
            IF (DQ0.GE.CSP) GOTO 160
            UP = UP - 1

! Scan column that matches with row Q0
            J = IPERM (Q0)
            VJ = DQ0 - A (JPERM (J) ) + U (Q0)
            DO K = IP (J), IP (J + 1) - 1
               I = IRN (K)
               IF (L (I) .GE.UP) CYCLE
! DNEW is new cost
               DNEW = VJ + A (K) - U (I)
! Do not update D(I) if DNEW ge cost of shortest path
               IF (DNEW.GE.CSP) CYCLE
               IF (IPERM (I) .EQ.0) THEN
! Row I is unmatched; update shortest path info
                  CSP = DNEW
                  ISP = K
                  JSP = J
               ELSE
! Row I is matched; do not update D(I) if DNEW is larger
                  DI = D (I)
                  IF (DI.LE.DNEW) CYCLE
                  IF (L (I) .GE.LOW) CYCLE
                  D (I) = DNEW
                  IF (DNEW.LE.DMIN) THEN
                     LPOS = L (I)
                     IF (LPOS.NE.0) CALL MC64FD (LPOS, QLEN, M, Q, D, L,&
                     2)
                     LOW = LOW - 1
                     Q (LOW) = I
                     L (I) = LOW
                  ELSE
                     IF (L (I) .EQ.0) THEN
                        QLEN = QLEN + 1
                        L (I) = QLEN
                     ENDIF
                     CALL MC64DD (I, M, Q, D, L, 2)
                  ENDIF
! Update tree
                  JJ = IPERM (I)
                  OUT (JJ) = K
                  PR (JJ) = J
               ENDIF
            END DO
         END DO

! If CSP = RINF, no augmenting path is found
  160    IF (CSP.EQ.RINF) GOTO 190
! Find augmenting path by tracing backward in PR; update IPERM,JPERM
         NUM = NUM + 1

!       write(14,*) 'NUM = ',NUM
!       write(14,*) 'Augmenting path found from unmatched col ',JORD

         I = IRN (ISP)
         IPERM (I) = JSP
         JPERM (JSP) = ISP
         J = JSP
         DO 170 JDUM = 1, NUM
            JJ = PR (J)
            IF (JJ.EQ. - 1) GOTO 180
            K = OUT (J)
            I = IRN (K)
            IPERM (I) = JJ
            JPERM (JJ) = K
            J = JJ
  170    END DO
! End of dummy loop; this point is never reached

! Update U for rows in Q(UP:N)
  180    DO 185 KK = UP, N
            I = Q (KK)
            U (I) = U (I) + D (I) - CSP
  185    END DO
  190    DO 191 KK = LOW, N
            I = Q (KK)
            D (I) = RINF
            L (I) = 0
  191    END DO
         DO 193 KK = 1, QLEN
            I = Q (KK)
            D (I) = RINF
            L (I) = 0
  193    END DO

      END DO
! End of main loop


! Set dual column variable in D(1:N)
 1000 DO 200 J = 1, N
         K = JPERM (J)
         IF (K.NE.0) THEN
            D (J) = A (K) - U (IRN (K) )
         ELSE
            D (J) = ZERO
         ENDIF
  200 END DO
      DO 1201 I = 1, M
         IF (IPERM (I) .EQ.0) U (I) = ZERO
 1201 END DO

      IF (NUM.EQ.N.AND.M.EQ.N) GOTO 1100
! Complete IPERM; L, JPERM are work arrays
      CALL MC64_HSL_X (M, N, IPERM, L, JPERM)

! The matrix is structurally singular, complete IPERM.
! JPERM, OUT are work arrays
!     DO 300 J = 1,N
!       JPERM(J) = 0
! 300 CONTINUE
!     K = 0
!     DO 310 I = 1,N
!       IF (IPERM(I).EQ.0) THEN
!         K = K + 1
!         OUT(K) = I
!       ELSE
!         J = IPERM(I)
!         JPERM(J) = I
!       ENDIF
! 310 CONTINUE
!     K = 0
!     DO 320 J = 1,N
!       IF (JPERM(J).NE.0) GO TO 320
!       K = K + 1
!       JDUM = OUT(K)
!       IPERM(JDUM) = - J
! 320 CONTINUE
 1100 RETURN
END SUBROUTINE MC64_HSL_W


!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_Z (M, N, IRN, LIRN, IP, LENC, IPERM, NUM, PR,&
      ARP, CV, OUT)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! PR(I) is the previous row to I in the depth first search.
!   It is used as a work array in the sorting algorithm.
! Elements (IPERM(I),I) I=1,...M  are non-zero at the end of the
!   algorithm unless N assignments have not been made.  In which case
!   (IPERM(I),I) will be zero for N-NUM entries.
! CV(I) is the most recent row extension at which column I was visited.
! ARP(I) is one less than the number of non-zeros in row I
!   which have not been scanned when looking for a cheap assignment.
! OUT(I) is one less than the number of non-zeros in row I
!   which have not been scanned during one pass through the main loop.
!
      INTEGER LIRN, M, N, NUM
      INTEGER ARP (N), CV (M), IRN (LIRN), IP (N), IPERM (M), LENC (N), &
      OUT (N), PR (N)

      INTEGER I, II, IN1, IN2, J, J1, JORD, K, KK

! Initialize variables to eliminate copmiler warnings
      II = -HUGE(II); IN2 = -HUGE(IN2)

      DO I = 1, M
         CV (I) = 0
         IPERM (I) = 0
      END DO
      DO J = 1, N
         ARP (J) = LENC (J) - 1
      END DO
      NUM = 0
!
! Main loop. Each pass round this loop either results in a new
! assignment or gives a row with no assignment.
!
      OUTER: DO JORD = 1, N
!
         J = JORD
         PR (J) = - 1
         SCAN: DO K = 1, JORD
! Look for a cheap assignment
            IN1 = ARP (J)
            IF (IN1.GE.0) THEN
               IN2 = IP (J) + LENC (J) - 1
               IN1 = IN2 - IN1
               DO II = IN1, IN2
                  I = IRN (II)
                  IF (IPERM (I) .EQ.0) EXIT scan
               END DO
! No cheap assignment in row.
               ARP (J) = - 1
            END IF
! Begin looking for assignment chain starting with row J.
            OUT (J) = LENC (J) - 1
! Inner loop.  Extends chain by one or backtracks.
            DO KK = 1, JORD
               IN1 = OUT (J)
               IF (IN1.GE.0) THEN
                  IN2 = IP (J) + LENC (J) - 1
                  IN1 = IN2 - IN1
! Forward scan.
                  DO II = IN1, IN2
                     I = IRN (II)
                     IF (CV (I) .EQ.JORD) CYCLE
! Column I has not yet been accessed during this pass.
                     J1 = J
                     J = IPERM (I)
                     CV (I) = JORD
                     PR (J) = J1
                     OUT (J1) = IN2 - II - 1
                     CYCLE SCAN
                  END DO
               ENDIF
! Backtracking step.
               J = PR (J)
               IF (J.EQ. - 1) CYCLE OUTER
            END DO
         END DO SCAN
!
! New assignment is made.
         IPERM (I) = J
         ARP (J) = IN2 - II - 1
         NUM = NUM + 1
         DO K = 1, JORD
            J = PR (J)
            IF (J.EQ. - 1) EXIT
            II = IP (J) + LENC (J) - OUT (J) - 2
            I = IRN (II)
            IPERM (I) = J
         END DO
!
      END DO OUTER


! IPERM is complete if M = N and NUM = N
      IF (M.EQ.N.and.NUM.EQ.N) RETURN

! Complete IPERM; CV, ARP are work arrays
      CALL MC64_HSL_X (M, N, IPERM, CV, ARP)

END SUBROUTINE MC64_HSL_Z

!**********************************************************************
!CCCC LAST UPDATE Tue Nov 26 03:20:26 MET 2002
SUBROUTINE MC64_HSL_X (M, N, IPERM, RW, CW)
      IMPLICIT NONE
!
! *** Copyright (c) 2002  I.S. Duff and J. Koster                   ***
! *** Although every effort has been made to ensure robustness and  ***
! *** reliability of the subroutines in this MC64 suite, we         ***
! *** disclaim any liability arising through the use or misuse of   ***
! *** any of the subroutines.                                       ***
!
! Complete the (incomplete) row permutation in IPERM.
!
      INTEGER M, N
      INTEGER RW (M), CW (N), IPERM (M)

      INTEGER I, J, K

! If M=N, the matrix is structurally singular, complete IPERM
! If M>N, the matrix is rectangular, complete IPERM

! RW, CW are work arrays;
! Store indices of unmatched rows in RW
! Mark matched columns in CW

      DO J = 1, N
         CW (J) = 0
      END DO
      K = 0
      DO I = 1, M
         IF (IPERM (I) .EQ.0) THEN
            K = K + 1
            RW (K) = I
         ELSE
            J = IPERM (I)
            CW (J) = I
         ENDIF
      END DO
      K = 0
      DO J = 1, N
         IF (CW (J) .NE.0) CYCLE
         K = K + 1
         I = RW (K)
         IPERM (I) = - J
      END DO
      DO J = N + 1, M
         K = K + 1
         I = RW (K)
         IPERM (I) = - J
      END DO
END SUBROUTINE MC64_HSL_X

end module hsl_mc64_double
! COPYRIGHT (c) 2007 Science & Technology Facilities Council
!
! Version: 3.2.0
! For version history see ChangeLog
!
    MODULE hsl_mc68_integer

      USE hsl_zb01_integer

      IMPLICIT NONE
      PRIVATE

! ---------------------------------------------------
! Precision
! ---------------------------------------------------
      INTEGER, PARAMETER :: myreal_mc68 = kind(1.0D0)
      INTEGER, PARAMETER :: myint = kind(1)
      INTEGER, PARAMETER :: long = selected_int_kind(18)

! ---------------------------------------------------
! Error flags
! ---------------------------------------------------
      INTEGER (myint), PARAMETER :: mc68_err_memory_alloc = -1, & ! memory
! alloc
! error
        mc68_err_memory_dealloc = -2, & ! memory dealloc error
        mc68_err_n = -3, & ! n<1
        mc68_err_ord = -4, & ! ord not associated with an ordering
        mc68_err_metis = -5, & ! MeTiS ordering requested but not linked
        mc68_err_zb01 = -6 ! None (de)allocation error from call to
! zb01_expand1

! ---------------------------------------------------
! Warning flags
! ---------------------------------------------------
      INTEGER (myint), PARAMETER :: mc68_warn_diag = 1, & ! No diags and ord=4
        mc68_warn_rank = 2, & ! Matrix rank deficient when using ord=4
        mc68_warn_rank_diag = 3 ! Matrix rank deficient, no diags and ord=4

! ---------------------------------------------------
! Derived type definitions
! ---------------------------------------------------
      TYPE, PUBLIC :: mc68_control
        INTEGER :: lp = 6 ! stream number for error messages
        INTEGER :: wp = 6 ! stream number for warning messages
        INTEGER :: mp = 6 ! stream number for diagnostic messages
        INTEGER :: nemin = 1 ! stream number for diagnostic messages
        INTEGER :: print_level = 0 ! amount of informational output required
        INTEGER :: row_full_thresh = 100 ! percentage threshold for full row
        INTEGER :: row_search = 10 ! Number of rows searched for pivot with
! ord=6
      END TYPE mc68_control

      TYPE, PUBLIC :: mc68_info
        INTEGER :: flag = 0 ! error/warning flag
        INTEGER :: iostat = 0 ! holds Fortran iostat parameter
        INTEGER :: stat = 0 ! holds Fortran stat parameter
        INTEGER :: out_range = 0 ! holds number of out of range entries
! ignored
        INTEGER :: duplicate = 0 ! holds number of duplicate entries
        INTEGER :: n_compressions = 0 ! holds number of compressions in order
        INTEGER :: n_zero_eigs = -1 ! holds the number of zero eigs from ma47
        INTEGER :: l_workspace = 0 ! holds length of workspace iw used in
! order
        INTEGER :: zb01_info = 0 ! holds flag from zb01_expand1 call
        INTEGER :: n_dense_rows = 0 ! holds number of dense rows from amdd
      END TYPE mc68_info

      INTERFACE mc68_order
        MODULE PROCEDURE mc68_order_integer
      END INTERFACE

      PUBLIC mc68_order

    CONTAINS

! ---------------------------------------------------------------

      SUBROUTINE mc68_order_integer(ord,n,ptr,row,perm,control,info, &
          min_l_workspace)
! subroutine mc68_order_integer(ord,A,perm,control,info) constructs
! an elimination order PERM for a symmetric matrix A using a chosen
! ordering ORD

! ord: is an INTEGER scalar with INTENT(in). It specifies which
! ordering is to be used. The choice is as follows:
! ord = 1 : Approximate minimum degree with provision for dense rows
! 2 : Minimum degree
! 3 : MeTiS
! 4 : MA47 ordering for indefinite matrices
        INTEGER, INTENT (IN) :: ord

! n: is an INTEGER scalar with INTENT(in). It must hold the number of rows in A
        INTEGER, INTENT (IN) :: n

! ptr: is an INTEGER array with INTENT(in) and size n. ptr(j) holds position in 
!       row of start of row indices for column j. ptr(n)+1 must equal the number
!       of entries stored + 1. Only the lower triangular entries are stored with
!       no duplicates or out-of-range entries
        INTEGER, INTENT (IN) :: ptr(n+1)

! row: is an INTEGER array with INTENT(in) and size at least as large as 
!       ptr(n+1)-1
        INTEGER, INTENT (IN) :: row(:)

! perm: is an integer array with intent(out) of size n. It holds
! the elimination order
        INTEGER (myint), INTENT (OUT) :: perm(n)

! control: is of derived type MC68_control with intent(in). Controls
! action
        TYPE (mc68_control), INTENT (IN) :: control

! info: is of derived type MC68_info with intent(out).
! info%flag
! = 0 if successful
! = MC68_ERR_MEMORY_ALLOC if memory allocation failed
! = MC68_ERR_MEMORY_DEALLOC if memory deallocation failed
! = MC68_ERR_N if n<1
! = MC68_ERR_METIS if MeTiS ordering is requested but not linked
        TYPE (mc68_info), INTENT (OUT) :: info

! min_l_workspace: is an optional integer scalare with intent(in). It
! specifies the minimum amount of workspace to be use
        INTEGER, INTENT (IN), OPTIONAL :: min_l_workspace

! ---------------------------------------------
! Local variables
! ---------------------------------------------
        INTEGER, ALLOCATABLE :: ipe(:) ! copy of pointers which is later
! modified
        INTEGER, ALLOCATABLE :: iw(:) ! copy of row indices
        INTEGER, ALLOCATABLE :: work1(:) ! work array
        INTEGER, ALLOCATABLE :: work2(:) ! work array
        INTEGER, ALLOCATABLE :: work3(:) ! work array
        INTEGER, ALLOCATABLE :: work4(:) ! work array
        INTEGER, ALLOCATABLE :: work6(:) ! work array
        INTEGER, ALLOCATABLE :: work7(:) ! work array
        INTEGER, ALLOCATABLE :: work8(:) ! work array
        INTEGER, ALLOCATABLE :: work9(:) ! work array
        INTEGER :: lp ! stream number for error messages
        INTEGER :: wp ! stream number for warning messages
        INTEGER :: mp ! stream number for diagnostic messages
        INTEGER (long) :: iwlen ! length of iw
        INTEGER :: iw1 ! work integers
        INTEGER :: i, k1, k2, iwfr, k, j, diag
        REAL (myreal_mc68) :: thresh
        LOGICAL :: printe ! errors to be printed?
        LOGICAL :: printw ! warnings to be printed?
        LOGICAL :: printi ! basic diagnostic to be printed?
        LOGICAL :: printd ! additional diagnostic to be printed?

! ---------------------------------------------
! Set stream numbers
! ---------------------------------------------
        lp = control%lp
        wp = control%wp
        mp = control%mp

! ---------------------------------------------
! Printing levels
! ---------------------------------------------
        printe = (control%print_level>=0 .AND. lp>=0)
        printw = (control%print_level>=0 .AND. wp>=0)
        printi = (control%print_level>=1 .AND. mp>=0)
        printd = (control%print_level>=2 .AND. mp>=0)
        IF (printi) THEN
          WRITE (mp,'(a)') ' '
          WRITE (mp,'(a)') 'MC68_order:'
        END IF

! Initialise info
        info%flag = 0
        info%stat = 0

! ---------------------------------------------
! Check that restrictions are adhered to
! ---------------------------------------------
        IF (n<1) THEN
          info%flag = mc68_err_n
          IF (printe) CALL mc68_print_message(info%flag,lp, &
            context='mc68_order')
          RETURN
        END IF

        IF (ord<1 .OR. ord>4) THEN
          info%flag = mc68_err_ord
          IF (printe) CALL mc68_print_message(info%flag,lp, &
            context='mc68_order')
          RETURN
        END IF

        SELECT CASE (ord)

        CASE (1)
          IF (printi) THEN
            WRITE (mp,'(a60)') 'Approximate minimum degree ordering'
          END IF

        CASE (2)
          IF (printi) THEN
            WRITE (mp,'(a60)') &
              'Minimum degree ordering using methodology of MA27'
          END IF

        CASE (3)
          IF (printi) THEN
            WRITE (mp,'(a60)') 'Nested bisection ordering using MeTiS'
          END IF

        CASE (4)
          IF (printi) THEN
            WRITE (mp,'(a60)') &
              'Ordering for indefinite matrices using methodology of MA47'
          END IF
        END SELECT

        IF (printi) THEN
          WRITE (mp,'(a,i15)') 'n  =  ', n
        END IF

        IF (n==1) THEN
! ---------------------------------------------
! Matrix is of order 1 so no need to call ordering subroutines
! ---------------------------------------------
          perm(1) = 1
          IF (printi) THEN
            WRITE (mp,'(a)') ' '
            WRITE (mp,'(a)') 'Matrix of order 1'
          END IF
          GO TO 20
        END IF

        IF (ord==3) THEN
! ---------------------------------------------
! Check whether MeTiS is linked
! ---------------------------------------------
          ALLOCATE (ipe(3),iw(2),work1(8),work2(2),STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_alloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF
          ipe = (/ 1, 2, 3 /)
          iw = (/ 1, 2 /)
          work1(1) = 0
          work1(2) = 3
          work1(3) = 1
          work1(4) = 2
          work1(5) = 0
          work1(6) = 1
          work1(7) = 200
          work1(8) = 1

          CALL metis_nodend(2,ipe,iw,1,work1,work2,perm(1:2))

          DEALLOCATE (ipe,iw,work1,work2,STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_dealloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

          IF (perm(1)==-1) THEN
            info%flag = mc68_err_metis
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF
        END IF

        SELECT CASE (ord)

        CASE (1)
! ----------------------------------
! Approximate minimum degree ordering with provision for
! dense rows
! ----------------------------------
          iw1 = 2*ptr(n+1)
          IF (present(min_l_workspace)) THEN
            iwlen = max(iw1+n,min_l_workspace)
          ELSE
            iwlen = iw1 + n
          END IF
          info%l_workspace = iwlen

          ALLOCATE (work1(n),work8(10),ipe(n+1),iw(iwlen),STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_alloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

! Fill iw and ipe removing any diagonal entries
          iw(:) = 0
          ipe(:) = 0

! Set ipe(j) to hold no. nonzeros in column j
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                ipe(i) = ipe(i) + 1
                ipe(j) = ipe(j) + 1
              END IF
            END DO
          END DO

! Set ipe(j) to point to where row indices will end in iw

          DO j = 2, n
            ipe(j) = ipe(j-1) + ipe(j)
          END DO
          ipe(n+1) = ipe(n) + 1

! Fill iw and ipe
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                iw(ipe(i)) = j
                iw(ipe(j)) = i
                ipe(i) = ipe(i) - 1
                ipe(j) = ipe(j) - 1
              END IF
            END DO
          END DO
          DO j = 1, n
            ipe(j) = ipe(j) + 1
          END DO
          iwfr = ipe(n+1) ! Index of next entry to be filled in iw

! Form ordering
          work1(1:n) = ipe(2:n+1) - ipe(1:n)
          work8(1) = 6
          work8(2) = 6
          work8(3) = -1
          work8(4) = 1 ! enable dense row detection
          work8(5) = huge(0)
          work8(6:10) = 0
          CALL amdd(n,iwlen,ipe,iwfr,work1,iw,perm,work8,info)

! Deallocate arrays
          DEALLOCATE (ipe,iw,work1,work8,STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_dealloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

        CASE (2)
! ----------------------------------
! Minimum degree ordering using methodology of MA27
! ----------------------------------

! Set length of iw
          iw1 = 2*ptr(n+1) - 1
          IF (present(min_l_workspace)) THEN
            iwlen = max(iw1+n,min_l_workspace)
          ELSE
            iwlen = iw1 + n
          END IF

! Allocate required arrays
          ALLOCATE (work2(n),ipe(n+1),iw(iwlen),STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_alloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

! Fill ipe and iw

! Set ipe(j) to hold no. nonzeros in column j
          ipe(:) = 0
          iw(:) = 0
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                ipe(i) = ipe(i) + 1
                ipe(j) = ipe(j) + 1
              END IF
            END DO
          END DO

! Set ipe(j) to point to where row indices will end in iw
          iw(1) = ipe(1)
          ipe(1) = ipe(1) + 1

          DO j = 2, n
            iw(ipe(j-1)+1) = ipe(j)
            ipe(j) = ipe(j-1) + ipe(j) + 1
          END DO
          ipe(n+1) = ipe(n) + 1

! Fill iw and ipe
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                iw(ipe(i)) = j
                iw(ipe(j)) = i
                ipe(i) = ipe(i) - 1
                ipe(j) = ipe(j) - 1
              END IF
            END DO
          END DO
          iwfr = ipe(n+1) ! Index of next entry to be filled in iw

          DO j = 1, n
            IF (iw(ipe(j))==0) ipe(j) = 0
          END DO



          thresh = float(control%row_full_thresh)/100.0
          iwfr = ipe(n+1)


          CALL mc68_min_deg_anal(n,ipe,iw,iwlen,iwfr,work2,huge(0), &
            info%n_compressions,thresh,info)


! set ipe correctly
          DO i = 1, n
            IF (work2(i)==0) THEN
              k1 = i
10            k2 = k1
              k1 = -ipe(k2)
              IF (work2(k1)==0) GO TO 10
              ipe(i) = -k1
            END IF
          END DO
          info%n_zero_eigs = -1

          CALL mc68_min_deg_tree_search(n,ipe,work2,perm,control%nemin,info)

          IF (info%flag<0 .AND. printe) THEN
            CALL mc68_print_message(info%flag,lp,context='mc68_order')
            RETURN
          END IF


! Deallocate required arrays
          DEALLOCATE (work2,ipe,iw,STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_dealloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

        CASE (3)
! ----------------------------------
! MeTiS ordering
! ----------------------------------

! Set length of iw
          IF (present(min_l_workspace)) THEN
            iwlen = max(2*ptr(n+1)-2,min_l_workspace)
          ELSE
            iwlen = 2*ptr(n+1) - 2
          END IF
          info%l_workspace = iwlen

! Allocate arrays
          ALLOCATE (work1(8),ipe(n+1),iw(iwlen),work3(n),STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_alloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

! Fill iw and ipe removing any diagonal entries
          iw(:) = 0
          ipe(:) = 0

! Set ipe(j) to hold no. nonzeros in column j
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                ipe(i) = ipe(i) + 1
                ipe(j) = ipe(j) + 1
              END IF
            END DO
          END DO

! Set ipe(j) to point to where row indices will end in iw

          DO j = 2, n
            ipe(j) = ipe(j-1) + ipe(j)
          END DO
          ipe(n+1) = ipe(n) + 1

! Fill iw and ipe
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                iw(ipe(i)) = j
                iw(ipe(j)) = i
                ipe(i) = ipe(i) - 1
                ipe(j) = ipe(j) - 1
              END IF
            END DO
          END DO
          DO j = 1, n
            ipe(j) = ipe(j) + 1
          END DO
          iwfr = ipe(n+1) ! Index of next entry to be filled in iw

! Carry out ordering
          work1(1) = 0
          work1(2) = 3
          work1(3) = 1
          work1(4) = 2
          work1(5) = 0
          work1(6) = 1
          work1(7) = 200
          work1(8) = 1
          CALL metis_nodend(n,ipe,iw,1,work1,work3,perm)

! Compression information not returned from meTiS
          info%n_compressions = 0
          info%n_zero_eigs = -1
          info%n_dense_rows = -1

! Deallocate arrays
          DEALLOCATE (ipe,iw,work1,work3,STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_dealloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

        CASE (4)

! ----------------------------------
! Ordering for indefinite matrices using methodology of MA47
! ----------------------------------
! Set length of iw
          iw1 = 2*ptr(n+1) - 2
          IF (present(min_l_workspace)) THEN
            iwlen = max(2*iw1+n,min_l_workspace)
          ELSE
            iwlen = 2*iw1 + n
          END IF

! Allocate required arrays
          ALLOCATE (work2(n),ipe(n+1),iw(iwlen),work3(n),work4(n),work6(n), &
            work7(n),work8(n),work9(n),STAT=info%stat)
          IF (info%stat/=0) THEN
            info%flag = mc68_err_memory_alloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

! Fill ipe and iw


! Set ipe(j) to hold no. nonzeros in column j
          ipe(:) = 0
          iw(:) = 0
          diag = 0
          work2 = -1
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                ipe(i) = ipe(i) + 1
                ipe(j) = ipe(j) + 1
              ELSE
                ipe(i) = ipe(i) + 1
                diag = diag + 1
                work2(i) = 1
              END IF
            END DO
          END DO

! Set ipe(j) to point to where row indices will end in iw
          iw(1) = ipe(1)
          ipe(1) = ipe(1) + 1

          DO j = 2, n
            iw(ipe(j-1)+1) = ipe(j)
            ipe(j) = ipe(j-1) + ipe(j) + 1
          END DO
          ipe(n+1) = ipe(n) + 1

! Fill iw and ipe
          DO j = 1, n
            DO k = ptr(j), ptr(j+1) - 1
              i = row(k)
              IF (j/=i) THEN
                iw(ipe(i)) = j
                iw(ipe(j)) = i
                ipe(i) = ipe(i) - 1
                ipe(j) = ipe(j) - 1
              ELSE
                iw(ipe(i)) = j
                ipe(i) = ipe(i) - 1
              END IF
            END DO
          END DO
          iwfr = ipe(n+1) ! Index of next entry to be filled in iw
! Warn if no diagonal entries present
          IF (diag==0) THEN
            info%flag = mc68_warn_diag
          END IF

! Carryout ordering
          CALL mc68_ma47_analyse(n,ipe,iw,iwlen,iwfr,control%row_search,work3, &
            work2,work4,info)
          info%l_workspace = iwlen
          info%n_dense_rows = -1
          IF (info%flag<0) THEN
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

          CALL mc68_ma47_treesearch(n,work4,work3,work2,perm,control%nemin, &
            info)

          IF (info%flag>0) THEN
            IF (printw) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
          END IF

          IF (info%flag<0) THEN
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF

! Deallocate required arrays
          DEALLOCATE (work3,work4,work6,work7,work8,work9,work2,ipe,iw, &
            STAT=info%stat)
          IF (info%stat<0) THEN
            info%flag = mc68_err_memory_dealloc
            IF (printe) CALL mc68_print_message(info%flag,lp, &
              context='mc68_order')
            RETURN
          END IF
        END SELECT


20      IF (printd) THEN
! ---------------------------------------------
! Print out perm
! ---------------------------------------------
          WRITE (mp,'(a7)') 'perm = '
          WRITE (mp,'(5i15)') (perm(i),i=1,n)
        ELSE IF (printi) THEN
! ---------------------------------------------
! Print out first few entries of perm
! ---------------------------------------------
          WRITE (mp,'(a21)') 'perm(1:min(5,n)) = '
          WRITE (mp,'(5i15)') (perm(i),i=1,min(5,n))
        END IF

        IF (printi) THEN
          CALL mc68_print_message(info%flag,mp,context='mc68_order')
        ELSE IF (printw .AND. (info%stat>0)) THEN
          CALL mc68_print_message(info%flag,wp,context='mc68_order')
        END IF

      END SUBROUTINE mc68_order_integer


      SUBROUTINE mc68_print_message(flag,unit,context)
! Prints out errors and warnings according to value of flag

! flag: is an integer scaler of intent(in). It is the information flag
! whose corresponding error message is printed
        INTEGER (myint), INTENT (IN) :: flag

! unit: is an integer scaler of intent(in). It is the unit number the
! error message should be printed on
        INTEGER (myint), INTENT (IN) :: unit

! context: is an optional assumed size character array of intent(in).
! It describes the context under which the error occured
        CHARACTER (len=*), OPTIONAL, INTENT (IN) :: context

        INTEGER (myint) :: length

        IF (unit<=0) RETURN

        IF (flag>0) THEN
          WRITE (unit,advance='yes',fmt='('' WARNING: '')')
        ELSE IF (flag<0) THEN
          WRITE (unit,advance='yes',fmt='('' ERROR: '')')
        END IF

        IF (present(context)) THEN
          length = len_trim(context)
          WRITE (unit,advance='no',fmt='('' '', a,'': '')') context(1:length)
        END IF

        SELECT CASE (flag)
        CASE (0)
          WRITE (unit,'(A)') 'successful completion'

        CASE (mc68_err_memory_alloc)
          WRITE (unit,'(A)') 'memory allocation failure'

        CASE (mc68_err_memory_dealloc)
          WRITE (unit,'(A)') 'memory deallocation failure'

        CASE (mc68_err_n)
          WRITE (unit,'(A)') 'restriction n>=1 violated'

        CASE (mc68_err_ord)
          WRITE (unit,'(A)') 'ord is not associated with an ordering'

        CASE (mc68_err_metis)
          WRITE (unit,'(A)') 'MeTiS ordering requested but not linked'

        CASE (mc68_err_zb01)
          WRITE (unit,'(A)') 'temporary file failure'

        CASE (mc68_warn_diag)
          WRITE (unit,'(A)') 'no diagonal entries'

        CASE (mc68_warn_rank)
          WRITE (unit,'(A)') 'matrix rank deficient'

        CASE (mc68_warn_rank+mc68_warn_diag)
          WRITE (unit,'(A)') 'no diagonal entries and matrix rank deficient'

        END SELECT

      END SUBROUTINE mc68_print_message




      SUBROUTINE mc68_min_deg_anal(n,ipe,iw,lw,iwfr,nv,iovflo,ncmpa,fratio, &
          info)
! F90 version of subroutine MA27HD.

! ANALYSIS SUBROUTINE

! GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
! PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
! IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
! VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
! I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
! TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
! OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).

! N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
! START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
! DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
! SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
! IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
! LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
! THE CREATED ELEMENT IS NULL. IF ELEMENT I
! IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
! IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
! ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
! DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
! LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
! ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
! IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
! IN THE NEW ELEMENT.
! LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
! IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
! NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
! JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
! THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
! VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
! IOVFLO should be set to a high legitimate integer.
! It is used as a flag.
! NCMPA number of compresses.
! FRATIO is the density of rows regarded as dense.

        REAL (myreal_mc68), INTENT (IN) :: fratio
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n, iovflo
        INTEGER, INTENT (INOUT) :: iwfr, ipe(n), iw(lw)
        INTEGER, INTENT (OUT) :: nv(n), ncmpa
        TYPE (mc68_info), INTENT (INOUT) :: info

! Local variables

! NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
! SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
! IF IT IS LAST IN ITS LIST.
! LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
! LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
! -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
! IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
! IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
! IF THERE ARE NONE.
! FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
! WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
! FLAG HAS THE FOLLOWING VALUES.
! A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
! FLAG(ME)=-1
! B) FOR VARIABLES JS
! FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
! FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
! FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
! ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
! CALCULATION
! FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
! ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
! CALCULATION
! C) FOR ELEMENTS IE
! FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
! FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
! CALCULATION FOR IS.
! FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
! DEGREE CALCULATION FOR IS
! ..
! .. Array Arguments ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: flag, ipd, lst, nxt

! ..
! .. Local Scalars ..
! LIMIT  Limit on number of variables for putting node in root.
! NVROOT Number of variables in the root node
! ROOT   Index of the root node (N+1 if none chosen yet).
        INTEGER :: i, id, idl, idn, ie, ip, is, jp, jp1, jp2, js, k, k1, k2, &
          ke, kp, kp0, kp1, kp2, ks, l, len, limit, ln, ls, lwfr, md, me, ml, &
          ms, nel, nflg, np, np0, ns, nvpiv, nvroot, root
! ..
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, min, nint
! ..
! If a column of the reduced matrix has relative density greater than
! CNTL(2), it is forced into the root. All such columns are taken to
! have sparsity pattern equal to their merged patterns, so the fill
! and operation counts may be overestimated.

! IS,JS,KS,LS,MS,NS ARE USED TO REFER TO SUPERVARIABLES.
! IE,JE,KE ARE USED TO REFER TO ELEMENTS.
! IP,JP,KP,K,NP ARE USED TO POINT TO LISTS OF ELEMENTS.
! OR SUPERVARIABLES.
! ID IS USED FOR THE DEGREE OF A SUPERVARIABLE.
! MD IS USED FOR THE CURRENT MINIMUM DEGREE.
! IDN IS USED FOR THE NO. OF VARIABLES IN A NEWLY CREATED ELEMENT
! NEL IS USED TO HOLD THE NO. OF VARIABLES THAT HAVE BEEN
! ELIMINATED.
! ME=MS IS THE NAME OF THE SUPERVARIABLE ELIMINATED AND
! OF THE ELEMENT CREATED IN THE MAIN LOOP.
! NFLG IS USED FOR THE CURRENT FLAG VALUE IN ARRAY FLAG. IT STARTS
! WITH THE VALUE IOVFLO AND IS REDUCED BY 1 EACH TIME IT IS USED
! UNTIL IT HAS THE VALUE 2 WHEN IT IS RESET TO THE VALUE IOVFLO.

! INITIALIZATIONS

        ALLOCATE (flag(n),ipd(n),lst(n),nxt(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF
        info%n_dense_rows = 0
        ipd(:) = 0
        nv(:) = 1
        flag(:) = iovflo
        md = 1
        ncmpa = 0
        nflg = iovflo
        nel = 0
        root = n + 1
        nvroot = 0

! LINK TOGETHER VARIABLES HAVING SAME DEGREE
        DO is = 1, n
          k = ipe(is)
          IF (k<=0) THEN
! WE HAVE A VARIABLE THAT CAN BE ELIMINATED AT ONCE BECAUSE THERE
! IS
! NO OFF-DIAGONAL NON-ZERO IN ITS ROW.
            nel = nel + 1
            flag(is) = -1
            nxt(is) = 0
            lst(is) = 0
          ELSE
            id = iw(k) + 1
            ns = ipd(id)
            IF (ns>0) lst(ns) = is
            nxt(is) = ns
            ipd(id) = is
            lst(is) = -id
          END IF
        END DO

! START OF MAIN LOOP
        DO ml = 1, n
! LEAVE LOOP IF ALL VARIABLES HAVE BEEN ELIMINATED.
          IF (nel+nvroot+1>=n) THEN
            GO TO 140
          ELSE
! FIND NEXT SUPERVARIABLE FOR ELIMINATION.
            DO id = md, n
              ms = ipd(id)
              IF (ms>0) GO TO 10
            END DO
10          md = id
! NVPIV HOLDS THE NUMBER OF VARIABLES IN THE PIVOT.
            nvpiv = nv(ms)

! REMOVE CHOSEN VARIABLE FROM LINKED LIST
            ns = nxt(ms)
            nxt(ms) = 0
            lst(ms) = 0
            IF (ns>0) lst(ns) = -id
            ipd(id) = ns
            me = ms
            nel = nel + nvpiv
! IDN HOLDS THE DEGREE OF THE NEW ELEMENT.
            idn = 0

! RUN THROUGH THE LIST OF THE PIVOTAL SUPERVARIABLE, SETTING TREE
! POINTERS AND CONSTRUCTING NEW LIST OF SUPERVARIABLES.
! KP IS A POINTER TO THE CURRENT POSITION IN THE OLD LIST.
            kp = ipe(me)
            flag(ms) = -1
! IP POINTS TO THE START OF THE NEW LIST.
            ip = iwfr
! LEN HOLDS THE LENGTH OF THE LIST ASSOCIATED WITH THE PIVOT.
            len = iw(kp)
            DO kp1 = 1, len
              kp = kp + 1
              ke = iw(kp)
! JUMP IF KE IS AN ELEMENT THAT HAS NOT BEEN MERGED INTO
! ANOTHER.
              IF (flag(ke)<=-2) THEN
! SEARCH VARIABLE LIST OF ELEMENT KE, USING JP AS A POINTER TO
! IT.
                ie = ke
                jp = ipe(ie)
                ln = iw(jp)
              ELSE
! JUMP IF KE IS AN ELEMENT THAT HAS BEEN MERGED INTO ANOTHER
! OR
! IS
! A SUPERVARIABLE THAT HAS BEEN ELIMINATED.
                IF (flag(ke)<=0) THEN
                  IF (ipe(ke)/=-root) THEN
                    GO TO 30
                  ELSE
! KE has been merged into the root
                    ke = root
                    IF (flag(ke)<=0) GO TO 30
                  END IF
                END IF
! WE HAVE A SUPERVARIABLE. PREPARE TO SEARCH REST OF LIST.
                jp = kp - 1
                ln = len - kp1 + 1
                ie = ms
              END IF

! SEARCH FOR DIFFERENT SUPERVARIABLES AND ADD THEM TO THE NEW
! LIST,
! COMPRESSING WHEN NECESSARY. THIS LOOP IS EXECUTED ONCE FOR
! EACH ELEMENT IN THE LIST AND ONCE FOR ALL THE SUPERVARIABLES
! IN THE LIST.
              DO jp1 = 1, ln
                jp = jp + 1
                is = iw(jp)
! JUMP IF IS IS NOT A PRINCIPAL VARIABLE OR HAS ALREADY BEEN
! COUNTED.
                IF (flag(is)<=0) THEN
                  IF (ipe(is)==-root) THEN
! IS has been merged into the root
                    is = root
                    iw(jp) = root
                    IF (flag(is)<=0) GO TO 20
                  ELSE
                    GO TO 20
                  END IF
                END IF
                flag(is) = 0
                IF (iwfr>=lw) THEN
! PREPARE FOR COMPRESSING IW BY ADJUSTING POINTERS AND
! LENGTHS SO THAT THE LISTS BEING SEARCHED IN THE INNER AND
! OUTER
! LOOPS CONTAIN ONLY THE REMAINING ENTRIES.
                  ipe(ms) = kp
                  iw(kp) = len - kp1
                  ipe(ie) = jp
                  iw(jp) = ln - jp1
! COMPRESS IW
                  CALL mc68_compress(n,ipe,iw,ip-1,lwfr,ncmpa)
! COPY NEW LIST FORWARD
                  jp2 = iwfr - 1
                  iwfr = lwfr
                  IF (ip<=jp2) THEN
                    DO jp = ip, jp2
                      iw(iwfr) = iw(jp)
                      iwfr = iwfr + 1
                    END DO
                  END IF
! ADJUST POINTERS FOR THE NEW LIST AND THE LISTS BEING
! SEARCHED.
                  ip = lwfr
                  jp = ipe(ie)
                  kp = ipe(me)
                END IF
! STORE IS IN NEW LIST.
                iw(iwfr) = is
                idn = idn + nv(is)
                iwfr = iwfr + 1
! REMOVE IS FROM DEGREE LINKED LIST
                ls = lst(is)
                lst(is) = 0
                ns = nxt(is)
                nxt(is) = 0
                IF (ns>0) lst(ns) = ls
                IF (ls<0) THEN
                  ls = -ls
                  ipd(ls) = ns
                ELSE IF (ls>0) THEN
                  nxt(ls) = ns
                END IF
20              CONTINUE
              END DO
! JUMP IF WE HAVE JUST BEEN SEARCHING THE VARIABLES AT THE END
! OF
! THE LIST OF THE PIVOT.
              IF (ie==ms) THEN
                GO TO 40
              ELSE
! SET TREE POINTER AND FLAG TO INDICATE ELEMENT IE IS ABSORBED
! INTO
! NEW ELEMENT ME.
                ipe(ie) = -me
                flag(ie) = -1
              END IF
30            CONTINUE
            END DO

! STORE THE DEGREE OF THE PIVOT.
40          nv(ms) = idn + nvpiv
! JUMP IF NEW ELEMENT IS NULL.
            IF (iwfr==ip) THEN
              ipe(me) = 0
            ELSE
              k1 = ip
              k2 = iwfr - 1

! RUN THROUGH NEW LIST OF SUPERVARIABLES REVISING EACH
! ASSOCIATED
! LIST,
! RECALCULATING DEGREES AND REMOVING DUPLICATES.
              limit = nint(fratio*(n-nel))
              DO k = k1, k2
                is = iw(k)
                IF (is/=root) THEN
                  IF (nflg<=2) THEN
! RESET FLAG VALUES TO +/-IOVFLO.
                    DO i = 1, n
                      IF (flag(i)>0) flag(i) = iovflo
                      IF (flag(i)<=-2) flag(i) = -iovflo
                    END DO
                    nflg = iovflo
                  END IF
! REDUCE NFLG BY ONE TO CATER FOR THIS SUPERVARIABLE.
                  nflg = nflg - 1
! BEGIN WITH THE DEGREE OF THE NEW ELEMENT. ITS VARIABLES
! MUST
! ALWAYS
! BE COUNTED DURING THE DEGREE CALCULATION AND THEY ARE
! ALREADY
! FLAGGED WITH THE VALUE 0.
                  id = idn
! RUN THROUGH THE LIST ASSOCIATED WITH SUPERVARIABLE IS
                  kp1 = ipe(is) + 1
! NP POINTS TO THE NEXT ENTRY IN THE REVISED LIST.
                  np = kp1
                  kp2 = iw(kp1-1) + kp1 - 1
                  DO kp = kp1, kp2
                    ke = iw(kp)
! TEST WHETHER KE IS AN ELEMENT, A REDUNDANT ENTRY OR A
! SUPERVARIABLE.
                    IF (flag(ke)==-1) THEN
                      IF (ipe(ke)/=-root) THEN
                        GO TO 60
                      ELSE
! KE has been merged into the root
                        ke = root
                        iw(kp) = root
                        IF (flag(ke)==-1) GO TO 60
                      END IF
                    END IF
                    IF (flag(ke)>=0) THEN
                      GO TO 70
                    ELSE
! SEARCH LIST OF ELEMENT KE, REVISING THE DEGREE WHEN
! NEW
! VARIABLES
! FOUND.
                      jp1 = ipe(ke) + 1
                      jp2 = iw(jp1-1) + jp1 - 1
                      idl = id
                      DO jp = jp1, jp2
                        js = iw(jp)
! JUMP IF JS HAS ALREADY BEEN COUNTED.
                        IF (flag(js)>nflg) THEN
                          id = id + nv(js)
                          flag(js) = nflg
                        END IF
                      END DO
! JUMP IF ONE OR MORE NEW SUPERVARIABLES WERE FOUND.
                      IF (id<=idl) THEN
! CHECK WHETHER EVERY VARIABLE OF ELEMENT KE IS IN NEW
! ELEMENT ME.
                        DO jp = jp1, jp2
                          js = iw(jp)
                          IF (flag(js)/=0) GO TO 50
                        END DO
! SET TREE POINTER AND FLAG TO INDICATE THAT ELEMENT
! KE
! IS ABSORBED
! INTO NEW ELEMENT ME.
                        ipe(ke) = -me
                        flag(ke) = -1
                        GO TO 60
                      END IF
! STORE ELEMENT KE IN THE REVISED LIST FOR SUPERVARIABLE
! IS AND FLAG IT.
50                    iw(np) = ke
                      flag(ke) = -nflg
                      np = np + 1
                    END IF
60                  CONTINUE
                  END DO
                  np0 = np
                  GO TO 90
! TREAT THE REST OF THE LIST ASSOCIATED WITH SUPERVARIABLE
! IS.
! IT
! CONSISTS ENTIRELY OF SUPERVARIABLES.
70                kp0 = kp
                  np0 = np
                  DO kp = kp0, kp2
                    ks = iw(kp)
                    IF (flag(ks)<=nflg) THEN
                      IF (ipe(ks)==-root) THEN
                        ks = root
                        iw(kp) = root
                        IF (flag(ks)<=nflg) GO TO 80
                      ELSE
                        GO TO 80
                      END IF
                    END IF
! ADD TO DEGREE, FLAG SUPERVARIABLE KS AND ADD IT TO NEW
! LIST.
                    id = id + nv(ks)
                    flag(ks) = nflg
                    iw(np) = ks
                    np = np + 1
80                  CONTINUE
                  END DO
! MOVE FIRST SUPERVARIABLE TO END OF LIST, MOVE FIRST
! ELEMENT
! TO END
! OF ELEMENT PART OF LIST AND ADD NEW ELEMENT TO FRONT OF
! LIST.
90                IF (id>=limit) THEN
! Treat IS as full. Merge it into the root node.
                    info%n_dense_rows = info%n_dense_rows + 1
                    IF (nvroot==0) THEN
                      root = is
                      ipe(is) = 0
                    ELSE
                      iw(k) = root
                      ipe(is) = -root
                      nv(root) = nv(root) + nv(is)
                      nv(is) = 0
                      flag(is) = -1
                    END IF
                    nvroot = nv(root)
                  ELSE
                    iw(np) = iw(np0)
                    iw(np0) = iw(kp1)
                    iw(kp1) = me
! STORE THE NEW LENGTH OF THE LIST.
                    iw(kp1-1) = np - kp1 + 1

! CHECK WHETHER ROW IS IS IDENTICAL TO ANOTHER BY LOOKING
! IN
! LINKED
! LIST OF SUPERVARIABLES WITH DEGREE ID AT THOSE WHOSE
! LISTS
! HAVE
! FIRST ENTRY ME. NOTE THAT THOSE CONTAINING ME COME FIRST
! SO THE
! SEARCH CAN BE TERMINATED WHEN A LIST NOT STARTING WITH
! ME
! IS
! FOUND.
                    js = ipd(id)
                    DO l = 1, n
                      IF (js<=0) THEN
                        GO TO 120
                      ELSE
                        kp1 = ipe(js) + 1
                        IF (iw(kp1)/=me) THEN
                          GO TO 120
                        ELSE
! JS HAS SAME DEGREE AND IS ACTIVE. CHECK IF
! IDENTICAL
! TO IS.
                          kp2 = kp1 - 1 + iw(kp1-1)
                          DO kp = kp1, kp2
                            ie = iw(kp)
! JUMP IF IE IS A SUPERVARIABLE OR AN ELEMENT NOT
! IN
! THE LIST OF IS.
                            IF (abs(flag(ie)+0)>nflg) GO TO 100
                          END DO
                          GO TO 110

100                       js = nxt(js)
                        END IF
                      END IF
                    END DO
! SUPERVARIABLE AMALGAMATION. ROW IS IS IDENTICAL TO ROW
! JS.
! REGARD ALL VARIABLES IN THE TWO SUPERVARIABLES AS BEING
! IN
! IS. SET
! TREE POINTER, FLAG AND NV ENTRIES.
110                 ipe(js) = -is
                    nv(is) = nv(is) + nv(js)
                    nv(js) = 0
                    flag(js) = -1
! REPLACE JS BY IS IN LINKED LIST.
                    ns = nxt(js)
                    ls = lst(js)
                    IF (ns>0) lst(ns) = is
                    IF (ls>0) nxt(ls) = is
                    lst(is) = ls
                    nxt(is) = ns
                    lst(js) = 0
                    nxt(js) = 0
                    IF (ipd(id)==js) ipd(id) = is
                    GO TO 130
! INSERT IS INTO LINKED LIST OF SUPERVARIABLES OF SAME
! DEGREE.
120                 ns = ipd(id)
                    IF (ns>0) lst(ns) = is
                    nxt(is) = ns
                    ipd(id) = is
                    lst(is) = -id
                    md = min(md,id)
                  END IF
                END IF
130             CONTINUE
              END DO

! RESET FLAGS FOR SUPERVARIABLES IN NEWLY CREATED ELEMENT AND
! REMOVE THOSE ABSORBED INTO OTHERS.
              DO k = k1, k2
                is = iw(k)
                IF (nv(is)/=0) THEN
                  flag(is) = nflg
                  iw(ip) = is
                  ip = ip + 1
                END IF
              END DO
              iwfr = k1
              flag(me) = -nflg
! MOVE FIRST ENTRY TO END TO MAKE ROOM FOR LENGTH.
              iw(ip) = iw(k1)
              iw(k1) = ip - k1
! SET POINTER FOR NEW ELEMENT AND RESET IWFR.
              ipe(me) = k1
              iwfr = ip + 1
            END IF
          END IF
        END DO

! Absorb any remaining variables into the root
140     DO is = 1, n
          IF (nxt(is)/=0 .OR. lst(is)/=0) THEN
            IF (nvroot==0) THEN
              root = is
              ipe(is) = 0
            ELSE
              ipe(is) = -root
            END IF
            nvroot = nvroot + nv(is)
            nv(is) = 0
          END IF
        END DO
! Link any remaining elements to the root
        DO ie = 1, n
          IF (ipe(ie)>0) ipe(ie) = -root
        END DO
        IF (nvroot>0) nv(root) = nvroot


        DEALLOCATE (flag,ipd,lst,nxt,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF


      END SUBROUTINE mc68_min_deg_anal


      SUBROUTINE mc68_min_deg_tree_search(n,ipe,nv,perm,nemin,info)

! Tree search
! Given son to father tree pointers, reorder so that eldest son has
! smallest degree and perform depth-first
! search to find pivot order and number of eliminations
! and assemblies at each stage.
! N must be set to the matrix order. It is not altered.
! IPE(I) must be set equal to -(father of node I) or zero if
! node is a root, if NV(I) > 0. If NV(I) = 0, then I is
! subordinate variable of a supervariable and -IPE(I) points to
! principal variable.  It is altered to point to its next
! younger brother if it has one, but otherwise is not changed.
! NV(I) must be set to zero if variable is a subordinate variable
! of a supervariable and to the degree otherwise.
! PERM is set to the new permutation after dfs of tree.  PERM(I) is
! the position of variable I in the pivot order.

        INTEGER, INTENT (IN) :: n, nemin
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        TYPE (mc68_info), INTENT (INOUT) :: info
        INTEGER, INTENT (OUT) :: perm(n)

! Local variables
! IPS(I) is used temporarily to hold
! -(eldest son of node I) if it has one and 0 otherwise. It is
! finally set to hold the position of node I in the order.
! NE(IS) is set to the number of variables
! eliminated at stage IS of the elimination.
! NA(IS) is set to the number of elements
! assembled at stage IS of the elimination.
! NODE (I) is used during the code
! to hold the number of subordinate variables for variable I and
! on output it holds
! the node (in dfs ordering) at which variable I is eliminated.
! It is also defined for subordinate variables.
! ND(IS) is set to the degree at stage IS of
! the elimination.
! NSTEPS is set to the number of elimination steps.
! NEMIN is used to control the amalgamation process between
! a son and its father (if the number of fully summed
! variables of both nodes is smaller than NEMIN).
! SUBORD(I) holds the first subordinate variable
! for variable I if I
! is a principal variable and holds the next subordinate variable
! if otherwise.  It is zero at the end of the chain.

! .. Array Arguments ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: fils, frere, ips, na, nd, ne, &
          node, subord
! ..
! .. Scalars ..
        INTEGER :: i, ib, if, ifson, il, in, inb, inf, infs, inl, ino, inos, &
          ins, insw, int, iperm, is, ison, k, l, nr, nr1, nsteps

! Initialisations

        ALLOCATE (fils(n),frere(n),ips(n),na(n),nd(n),ne(n),node(n),subord(n), &
          STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

        ips(:) = 0
        ne(:) = 0
        node(:) = 0
        subord(:) = 0


! Set IPS(I) to -(eldest son of node I) and IPE(I) to next younger
! brother of node I if it has one.
        nr = n + 1
        DO i = 1, n
          if = -ipe(i)
          IF (nv(i)==0) THEN
! I is a subordinate node, principal variable is IF
            IF (subord(if)/=0) subord(i) = subord(if)
            subord(if) = i
            node(if) = node(if) + 1
! Node IF is the father of node I.
          ELSE IF (if/=0) THEN
! IS is younger brother of node I.
! IPS(IF) will eventually point to - eldest son of IF.
            is = -ips(if)
            IF (is>0) ipe(i) = is
            ips(if) = -i
          ELSE
! I is a root node
            nr = nr - 1
            ne(nr) = i
          END IF
        END DO

! We reorganize the tree so that the eldest son has maximum number of
! variables.  We combine nodes when the number of variables in a son
! is greater than or equal to the number of variables in the father.
! If the eldest son has the maximum number of variables,
! and if a combination is possible, it has to be possible with
! the eldest son.
! FILS is just used as workspace during this reorganization and is
! reset
! afterwards.

        DO i = 1, n
          fils(i) = ips(i)
        END DO
        nr1 = nr
        ins = 0
10      CONTINUE
! Jump if all roots processed.
        IF (nr1<=n) THEN
! Get next root
          ins = ne(nr1)
          nr1 = nr1 + 1
20        CONTINUE
! Depth first search through eldest sons.
          inl = fils(ins)
          IF (inl<0) THEN
            ins = -inl
            GO TO 20
          ELSE
            DO WHILE (ipe(ins)<0)
! INS is youngest son otherwise IPE value would be positive.
              ins = -ipe(ins)
! INS is now the father of the reorganized son so we can
! clear the pointer to the sons.
              fils(ins) = 0
! Continue backtracking until we encounter node with younger
! brother.
            END DO

            IF (ipe(ins)/=0) THEN
! INB is younger brother of INS.
              inb = ipe(ins)
              IF (nv(inb)>=nv(ins)) THEN
                ins = inb
! Do depth first search from younger brother
              ELSE
! Exchange INB and INS
! Find previous brother of INS (could be the father)
! then we do depth first search with INS = INB
                inf = inb
30              CONTINUE
                inf = ipe(inf)
                IF (inf>0) GO TO 30
! -INF IS THE FATHER
                inf = -inf
                infs = -fils(inf)
! INFS is eldest son of INF
                IF (infs==ins) THEN
! INS is eldest brother .. a role which INB now assumes
                  fils(inf) = -inb
                  ips(inf) = -inb
                  ipe(ins) = ipe(inb)
                  ipe(inb) = ins
                ELSE
                  insw = infs
40                CONTINUE
                  infs = ipe(insw)
                  IF (infs/=ins) THEN
                    insw = infs
                    GO TO 40
                  END IF
                  ipe(ins) = ipe(inb)
                  ipe(inb) = ins
                  ipe(insw) = inb
                END IF
                ins = inb
! Depth first search from moved younger brother
              END IF
              GO TO 20
            END IF
          END IF
! INS is a root, check for next one.
          ins = 0
          GO TO 10
        END IF
! Set FRERE and FILS
        DO i = 1, n
          frere(i) = ipe(i)
          fils(i) = ips(i)
        END DO

! Depth-first search.
! IL holds the current tree level. Roots are at level N, their sons
! are at level N-1, etc.
! IS holds the current elimination stage. We accumulate the number
! of eliminations at stage is directly in NE(IS). The number of
! assemblies is accumulated temporarily in NA(IL), for tree
! level IL, and is transferred to NA(IS) when we reach the
! appropriate stage IS.
        is = 1
! I is the current node.
        i = 0
! IPERM is used as pointer to setting permutation vector
        iperm = 1
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
! Stop if all roots used (needed because of subordinate variables)
            IF (nr>n) THEN
              GO TO 130
            ELSE
              i = ne(nr)
              ne(nr) = 0
              nr = nr + 1
              il = n
              na(n) = 0
            END IF
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in IPS as each is used and setting NA(IL)=0 for all levels
! reached.
          DO l = 1, n
            IF (ips(i)>=0) THEN
              GO TO 50
            ELSE
              ison = -ips(i)
              ips(i) = 0
              i = ison
              il = il - 1
              na(il) = 0
            END IF
          END DO
! Record position of node I in the order.
50        ips(i) = k
! Add number of subordinate variables to variable I
          ne(is) = ne(is) + node(i) + 1
          IF (il<n) na(il+1) = na(il+1) + 1
          na(is) = na(il)
          nd(is) = nv(i)
          node(i) = is
          perm(i) = iperm
          iperm = iperm + 1
! Order subordinate variables to node I
          in = i
60        CONTINUE
          IF (subord(in)/=0) THEN
            in = subord(in)
            node(in) = is
            perm(in) = iperm
            iperm = iperm + 1
            GO TO 60
          END IF
! Check for static condensation
          IF (na(is)==1) THEN
            IF (nd(is-1)-ne(is-1)==nd(is)) GO TO 70
          END IF
! Check for small numbers of eliminations in both last two steps.
          IF (ne(is)<nemin) THEN
            IF (na(is)/=0) THEN
              IF (ne(is-1)<nemin) GO TO 70
            END IF
          END IF
          is = is + 1
          GO TO 120

! Combine the last two steps
70        na(is-1) = na(is-1) + na(is) - 1
          nd(is-1) = nd(is) + ne(is-1)
          ne(is-1) = ne(is) + ne(is-1)
          ne(is) = 0
          node(i) = is - 1
! Find eldest son (IFSON) of node I (IS)
! Note that node I must have a son (node IS-1 is youngest)
          ifson = -fils(i)
! Now find youngest son INO (he is node IS-1)
          in = ifson
80        CONTINUE
          ino = in
          in = frere(in)
          IF (in>0) GO TO 80
! Cannot be root node .. so points to father
! Merge node IS-1 (INO) into node IS (I)
          nv(ino) = 0
! IPE already set .. was father pointer now principal variable
! pointer
! Now make subsidiary nodes of INO into subsidiary nodes of I.
! Subordinate nodes of INO become subordinate nodes of I
          in = i
90        CONTINUE
          IF (subord(in)/=0) THEN
            in = subord(in)
            node(in) = is - 1
            GO TO 90
          END IF
          subord(in) = ino
          in = ino
          IF (subord(in)/=0) THEN
            in = subord(in)
            ipe(in) = -i
          END IF

! INOS is eldest son of INO
          inos = -fils(ino)

! Find elder brother of node INO
! First check to see if he is only son
          IF (ifson/=ino) THEN
            in = ifson
100         CONTINUE
            ins = in
            in = frere(in)
            IF (in/=ino) GO TO 100
! INS is older brother .. make him brother of first son of INO (ie
! INOS)
! and  make INOS point to I now as father.
! Jump if there is no son of INO
            IF (inos==0) THEN
! Elder brother of INO just points to (new) father.
              frere(ins) = -i
              GO TO 120
            ELSE
              frere(ins) = inos
            END IF
          END IF
! INT is youngest brother of INOS.  Make him point to (new) father.
          in = inos
          IF (in/=0) THEN
110         CONTINUE
            int = in
            in = frere(in)
            IF (in>0) GO TO 110
            frere(int) = -i
          END IF
120       ib = ipe(i)
          IF (ib>=0) THEN
! Node I has a younger brother or is a root
            IF (ib>0) na(il) = 0
            i = ib
          ELSE
! I has no brothers. Go to father of node I
            i = -ib
            il = il + 1
          END IF
        END DO
130     nsteps = is - 1


        DEALLOCATE (fils,frere,ips,na,nd,ne,node,subord,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

      END SUBROUTINE mc68_min_deg_tree_search

      SUBROUTINE mc68_compress(n,ipe,iw,lw,iwfr,ncmpa)
! Is identical to subroutine MA27UD.  Internal version for MA57.
! COMPRESS LISTS HELD BY MA27H/HD (MA57H/HD) IN IW AND ADJUST POINTERS
! IN IPE TO CORRESPOND.
! N IS THE MATRIX ORDER. IT IS NOT ALTERED.
! IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
! ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
! IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
! LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
! LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
! IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
! LOCATION IN IW.
! ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
! NCMPA is number of compresses.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lw, n
        INTEGER, INTENT (OUT) :: iwfr
        INTEGER, INTENT (INOUT) :: ncmpa
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: ipe(n), iw(lw)
! ..
! .. Local Scalars ..
        INTEGER i, ir, k, k1, k2, lwfr
! ..
        ncmpa = ncmpa + 1
! PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
! LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
! -(LIST NUMBER).
        DO i = 1, n
          k1 = ipe(i)
          IF (k1>0) THEN
            ipe(i) = iw(k1)
            iw(k1) = -i
          END IF
        END DO

! COMPRESS
! IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
! LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
        iwfr = 1
        lwfr = iwfr
        DO ir = 1, n
          IF (lwfr>lw) THEN
            GO TO 20
          ELSE
! SEARCH FOR THE NEXT NEGATIVE ENTRY.
            DO k = lwfr, lw
              IF (iw(k)<0) GO TO 10
            END DO
            GO TO 20
! PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW
! POINTER
! AND PREPARE TO COPY LIST.
10          i = -iw(k)
            iw(iwfr) = ipe(i)
            ipe(i) = iwfr
            k1 = k + 1
            k2 = k + iw(iwfr)
            iwfr = iwfr + 1
            IF (k1<=k2) THEN
! COPY LIST TO NEW POSITION.
              DO k = k1, k2
                iw(iwfr) = iw(k)
                iwfr = iwfr + 1
              END DO
            END IF
            lwfr = k2 + 1
          END IF
        END DO
20      RETURN

      END SUBROUTINE mc68_compress
! -----------------------------------------------------------

      SUBROUTINE mc68_ma47_analyse(n,ipe,iw,lw,iwfr,cntl,count,nv,next,info)
        INTEGER, PARAMETER :: wp = kind(0.0D0)

! Analysis subroutine

! Given representation of the whole matrix, perform Markowitz
! ordering,
! constructing tree pointers.  It works with supervariables which
! are collections of one or more variables whose rows have the same
! pattern, starting with supervariable I containing variable I for I
! = 1,2,...,N.  Each supervariable has as numerical name that of one
! of its variables (its principal variable). When a supervariable is
! eliminated, its name is used for the name of the element created.


! N must be set to the matrix order. It is not altered.
! Constraint: N <= HUGE/3.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N.
! During execution, IPE(I) is used as follows:
! if I is a supervariable, IPE(I) is the position in IW of its list;
! if I is an element, IPE(I) is the position in IW of its list, or
! 0 if the list has zero length;
! if I is a variable that has been absorbed in a supervariable,
! IPE(I) = 0.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. The lists include any diagonal entries
! and the entries of both the upper and lower triangular parts.
! During execution, it holds lists for supervariables and lists for
! elements. A list for a supervariable represents its row of the
! reduced matrix and has the form:
! length ( = total number of elements and supervariables)
! list of element parts in the form JS+N*PART, where JS is an
! element name and PART has the value
! 2 for the part with leading zeros;
! 0 for the full part;
! 1 for the part with trailing zeros;
! list of supervariables.
! A list for an element has the form:
! length ( = 2 + total number of supervariables)
! number of supervariables in the part with leading zeros
! number of supervariables in the part with trailing zeros
! list of supervariables.
! LW must be set to the length of IW. It may be altered.
! IWFR must be set to the position in IW of the first unused location.
! It is revised during execution and continues to have
! this meaning.
! COUNT need not be set on entry. For supervariables, it is used to
! hold the row counts in the reduced matrix. COUNT(I) is negated
! if supervariable I has been eliminated as the trailing
! part of a tile or oxo pivot.
! NV must be set on entry so that NV(I) = 1 for a nondefective
! variable
! (nonzero diagonal entry) and NV(I) = -1 otherwise, I = 1,..., N.
! During execution, NV(JS) is used as follows:
! if JS is a variable that has been absorbed into another,
! NV(JS) = 0;
! if JS is a supervariable, NV(JS) holds its number of
! variables, negated for a defective supervariable.
! NEXT need not be set. During execution, it is used as follows:
! if JS is a supervariable that has not been eliminated or
! absorbed, NEXT(JS) is the next such supervariable having the
! same row count, or zero if JS is last in its list;
! if supervariable IS has been absorbed and JS is the next
! variable in its supervariable, NEXT(IS)=-JS;
! if supervariable IS was eliminated in a block pivot with
! supervariable JS, NEXT(IS)=-JS;
! if IE is an element whose parent in the tree is JE, NEXT(IE)=-JE;
! if IE is an element with no parent in the tree, NEXT(IE)=0.
! NCMP is number of times garbage collection routine is called.
! ICNTL must be set by the user as follows and is not altered.
! ICNTL(1)  must be set to the stream number for error messages.
! A value less than 1 suppresses output.
! ICNTL(2) must be set to the stream number for diagnostic output.
! A value less than 1 suppresses output.
! ICNTL(3) must be set to control the amount of output:
! 0 None.
! 1 Error messages only.
! 2 Error and warning messages.
! 3 As 2, plus scalar parameters and a few entries of array
! parameters on entry and exit.
! 4  As 2, plus all parameters on entry and exit.
! CNTL must set to one if the user wants the pivot order chosen
! by Markowitz strategy,, and greater than 1 for Markowits strategy
! with
! each search for a structured pivot limited to this number of rows.
! INFO is left unaltered except
! 8 is added to INFO%flag and
! INFO%n_zero_eigs is set to the number of zero eigenvalues found.

! Local constants

! Local variables
! BOUND True if row count is to be bounded rather than recomputed.
! CMIN holds the minimum Markowitz cost of a potential structured
! pivot
! found so far.
! COST holds the Markowitz cost of a potential structured pivot.
! FLG Flag value.
! FLG1 Flag value.
! I Temporary variable.
! IE is an element.
! IEL is an element.
! IP is used to point to a list of elements or supervariables.
! IR is used for row counts.
! IRL last value of IR
! IS is a supervariable.
! ITHR Index for do loop in case threshold needs to be raised.
! JP is used to point to a list of elements or supervariables.
! JP1 start of range for JP.
! JP2 end of range for JP.
! JP3 start of range for JP.
! JP4 end of range for JP.
! JS is a supervariable.
! K is used to point to a list of elements or supervariables.
! KE is an element.
! KIND has the value
! 1 for a full pivot,
! 2 for a tile pivot, or
! 3 for an oxo pivot
! KP is used to point to a list of elements or supervariables.
! KP1 start of range for KP.
! KP2 end of range for KP.
! KS is a supervariable.
! K1 start of range for K.
! K2 end of range for K.
! LIST is used to search a list.
! LOOP Row of the pivot when constructing the index list of a new
! element.
! LS is a supervariable.
! ME is the element created by the current pivot step.
! MINF is lower bound for the row count of a nondefective variable.
! MINR is used for the current minimum row count.
! ML is used for the main loop index.
! MP is stream for warning and diagnostic messages.
! MROWS holds the max. number of rows to be searched for a structured
! pivot.
! MS is a supervariable.
! NEL is used to hold the number of variables that have been
! eliminated.
! NFLG is used for the current flag value in array FLAG. It starts
! with the value N*3 and is reduced by 1 each time it is used
! until it has the value 4 when it is reset to the value N*3.
! NP is used to point to a list of elements or supervariables.
! NP0 is used to point to a list of elements or supervariables.
! NR is used for the row counts in the generated element:
! NR(1) for a row of the middle block
! NR(2) for a row of the last block
! NR(3) for a row of the first block
! NROWS holds the number of rows searched for a structured pivot.
! NS is a supervariable.
! NSC is used for supervariable counts:
! NSC(1) first zero block
! NSC(2) full block
! NSC(3) second zero block
! NSVARS Supervariable counts.
! NVC Number of variables in both pivot rows, excluding the pivot.
! NVPIV holds the number of variables in the pivot (each half in the
! case of a tile or oxo).
! PART Part of element: 0 for full part, 1 for trailing part,
! 2 for leading part.
! PIV holds the name(s) of the pivot supervariable(s), ordered in the
! case of a tile pivot so that the first diagonal entry is nonzero.
! PIVOT As PIV, except that this holds the original names when one
! has been split.
! PIVT In the case of a structured pivot with unbalanced
! supervariables,
! PIVT is equal to whichever of PIVOT(1) and PIVOT(2) has more
! variables. Zero otherwise.
! RANK Upper bound for the rank.
! SHIFT Movement of current list after compress.
! THRESH Row counts greater than THRESH are lower bounds (fill-ins may
! have been ignored.
! LAST need not be set. During execution, it is used as follows:
! if JS is a supervariable that has not been eliminated or
! absorbed, LAST(JS) is the previous such supervariable having
! the same row count, or zero if JS is first in its list;
! otherwise, if NEXT(JS) < 0, LAST(JS) is one of the values
! NEXT(JS), NEXT(-NEXT(JS)), NEXT(-NEXT(-NEXT(JS))), ...
! IPR need not be set. During execution IPR(IR) is the first
! supervariable with row count IR or zero if there are none.
! FLAG is used as workspace for element and supervariable flags.
! if IS is a pivotal supervariable, FLAG(IS) = -1;
! if IS is a non-pivotal supervariable, FLAG(IS) >= 0;
! if IS is a supervariable involved in the current pivot row (other
! than in the pivot):
! for a full pivot:           FLAG(IS) = 1
! for a tile or oxo pivot:
! if in both rows,       FLAG(IS) = 0
! if in first row only,  FLAG(IS) = 1
! if in second row only, FLAG(IS) = 2
! if IE is an element, FLAG(IE) < -1;
! if I is neither a supervariable nor an element, FLAG(I) = -1.
! LEAF need not be set. During execution, if IS is a supervariable,
! its
! final variable is LEAF(IS). LEAF(IS) is linked to IS through NEXT
! in a chain that covers the other variables of the supervariable.
! For an element IE to which the current supervariable IS belongs,
! LEAF(IE) holds the part to which the supervariable belongs.
! SVARS is used as workspace for the supervariables of the new
! element.
! .. Parameters ..
        REAL (wp) :: zero
        PARAMETER (zero=0E0_wp)
! ..
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n, cntl
        INTEGER, INTENT (INOUT) :: iwfr
        INTEGER (long), INTENT (INOUT) :: lw

! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        INTEGER, DIMENSION (:), ALLOCATABLE, INTENT (INOUT) :: iw
        INTEGER, INTENT (OUT) :: count(n), next(n)

        TYPE (mc68_info), INTENT (INOUT) :: info
! ..
! .. Local Scalars ..

        REAL (wp) :: cmin, cost

        INTEGER :: flg, flg1, ie, iel, ip, ir, irl, is, ithr, jp, jp1, jp2, &
          jp3, jp4, js, k, k1, k2, ke, kind1, kp, kp1, kp2, ks, list, loop, &
          ls, me, minf, minr, ml, mp, mrows, ms, nel, nflg, np, np0, nrows, &
          ns, nsvars, nvc, nvpiv, part, pivt, rank, shift, thresh
        INTEGER (long) :: szw
        LOGICAL :: bound
! ..
! .. Local Arrays ..
        INTEGER :: nr(3), nsc(3), piv(2), pivot(2), icntl4
        INTEGER, ALLOCATABLE, DIMENSION (:) :: flag, ipr, last, leaf, svars

        TYPE (zb01_info) :: info_zb01
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, max, min, real, sign, sqrt

        IF (cntl<=1) THEN
          icntl4 = 0
        ELSE
          icntl4 = cntl
        END IF
        mp = 6
        mrows = n
        IF (icntl4>1) mrows = icntl4

        ALLOCATE (flag(n),ipr(n),last(n),leaf(n),svars(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

! Create supervariables from variables with identical rows
        CALL mc68_ma47_merge(n,ipe,iw,lw,nv,next,last,leaf,flag,count,svars)

! Initializations

        thresh = max(sqrt(real(n,kind=wp)),3._wp)
        rank = n
        minf = 1
        minr = 1
        info%n_compressions = 0
        nflg = n*3
        nel = 0
        DO is = 1, n
          ipr(is) = 0
          flag(is) = nflg
        END DO

! Link together supervariables having same row count
        DO is = 1, n
          k = ipe(is)
          IF (k==0) THEN
! Row IS is identical to another row and has been absorbed.
            flag(is) = -1
          ELSE
            ir = iw(k)
            IF (ir==0) THEN
              rank = rank + nv(is)
! Treat zero block row as if it had nonzero diagonal entries.
              nv(is) = -nv(is)
              ir = 1
            END IF
! Store the row count and add the row to its linked list.
            count(is) = ir
            ns = ipr(ir)
            IF (ns>0) last(ns) = is
            next(is) = ns
            ipr(ir) = is
            last(is) = 0
          END IF
        END DO

! Start of main loop
        DO ml = 1, n
! Leave loop if all variables have been eliminated.
          IF (nel>=n) THEN
            GO TO 280
          ELSE
! Find minimum row count
            ir = minr
            DO minr = ir, n
              IF (ipr(minr)/=0) GO TO 10
            END DO
! Find next pivot.
10          DO ithr = 1, n
! Outer loop on ITHR needed in case of the threshold needing to
! be
! increased.
              nrows = 0
              cmin = real(n,kind=wp)**2
              DO ir = minr, n
                IF (minf<=ir) THEN
! Look for a full pivot
                  ms = ipr(ir)
                  DO list = 1, n
                    IF (ms==0) THEN
                      GO TO 20
! If this is a full pivot, accept it or raise the
! threshold and try
! again.
                    ELSE IF (nv(ms)>0) THEN
                      GO TO 60
                    ELSE
                      ms = next(ms)
                    END IF
                  END DO
                  GO TO 30
20                minf = ir + 1
                END IF
! Look for a tile or oxo pivot with MS as its defective
! variable.
30              ms = ipr(ir)
                DO list = 1, n
                  IF (ms==0) THEN
                    GO TO 50
                  ELSE
! We need not look for a tile with MS as its non-defective
! variable
! because if the row count of its mate is IS, its M-cost
! is
! (IR+IS-3)*(IS-1)
! >= (IS-1)**2 since IR >= 2
! >= (IR-1)**2 since IS >= IR
! Reduce NFLG by one to cater for this supervariable.
                    IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
                    nrows = nrows + 1
                    nflg = nflg - 1
                    kp = ipe(ms)
                    kp1 = kp + 1
                    kp2 = kp + iw(kp)
                    DO kp = kp1, kp2
                      part = (iw(kp)-1)/n
                      ke = iw(kp) - part*n
                      IF (flag(ke)/=-1) THEN
                        IF (flag(ke)<=-2) THEN
! KE is an element touched by MS.
                          jp = ipe(ke)
                          IF (part==2) THEN
! Supervariable in leading part of element
                            jp1 = jp + 3 + iw(jp+1)
                            jp2 = jp + iw(jp)
                          ELSE
! Supervariable in trailing part of element
                            jp1 = jp + 3
                            jp2 = jp + iw(jp) - iw(jp+2)
                          END IF
                        ELSE
! Supervariable MS is defective so cannot lie in
! full
! part of element
! We have reached the list of variables
                          jp1 = kp
                          jp2 = kp2
                        END IF
! Search the variable list.
                        DO jp = jp1, jp2
                          is = iw(jp)
                          IF (flag(is)>nflg) THEN
! This is the first time IS has been encountered
! in
! this search.
                            flag(is) = nflg
                            IF (nv(is)<0) THEN
! Potential oxo pivot
                              cost = real(count(is)-1,kind=wp)* &
                                real(ir-1,kind=wp)
                            ELSE
! Potential tile pivot
                              cost = real(ir+count(is)-3,kind=wp)* &
                                real(ir-1,kind=wp)
                            END IF

                            IF (cost<cmin) THEN
! This is the best potential pivot so far found.
                              cmin = cost
                              pivot(1) = is
                              pivot(2) = ms
! Exit loop no better pivot is possible.
                              IF (cmin<=real(ir-1,kind=wp)**2) GO TO 70
                            END IF
                          END IF
                        END DO

! Exit loop if we have searched the list of variables.
                        IF (jp2==kp2) GO TO 40
                      END IF
                    END DO
! Exit loop if enough rows have been searched.
40                  IF (nrows>=mrows) THEN
                      GO TO 70
                    ELSE
                      ms = next(ms)
                    END IF
                  END IF
                END DO
! Exit loop no better pivot is possible.
50              IF (cmin<=real(ir,kind=wp)**2) GO TO 70
              END DO
              GO TO 70
60            IF (ir<=thresh) THEN
                GO TO 90
              ELSE
                GO TO 80
              END IF
70            ir = max(count(pivot(1)),count(pivot(2)))
              IF (ir<=thresh) GO TO 100
! Revise the threshold and repeat the choice of pivot.
80            CALL mc68_ma47_thresh(thresh,ir+n/10,n,ipe,iw,lw,count,nv,next, &
                last,ipr,flag,nflg)
            END DO

! Accept full pivot
90          kind1 = 1
            me = ms
            pivot(1) = ms
            pivot(2) = ms
            pivt = 0
            nvpiv = nv(ms)
            cmin = real(ir-1,kind=wp)**2
            flag(ms) = -1
            nel = nel + nvpiv
            GO TO 110

! Accept tile or oxo pivot
! PIVOT(2) cannot have greater row count than PIVOT(1).
100         kind1 = 2
            IF (nv(pivot(1))<0) kind1 = 3
            flag(pivot(1)) = -1
            flag(pivot(2)) = -1
            piv(1) = pivot(1)
            piv(2) = pivot(2)
            nvpiv = abs(nv(pivot(1)))
            IF (nvpiv==abs(nv(pivot(2)))) THEN
              pivt = 0
            ELSE
              IF (nvpiv>abs(nv(pivot(2)))) THEN
                pivt = pivot(1)
                nvpiv = abs(nv(pivot(2)))
              ELSE
                pivt = pivot(2)
              END IF
! Split PIVT
              ms = leaf(pivt)
              DO k = 2, nvpiv
                ms = -next(ms)
              END DO
              leaf(ms) = leaf(pivt)
              leaf(pivt) = -next(ms)
! PIVT remains as an active supervariable, but with less
! variables.
              flag(pivt) = nflg
              next(ms) = 0
              nv(pivt) = sign(abs(nv(pivt))-nvpiv,nv(pivt))
              nv(ms) = sign(nvpiv,nv(pivt))
              count(ms) = count(pivt)
              IF (pivt==pivot(1)) THEN
                piv(1) = ms
              ELSE
                piv(2) = ms
              END IF
            END IF

            me = piv(1)
            nsc(2) = 0
            nvc = 0
            nel = nel + nvpiv*2

! Find variables of new element
110         nsvars = 0
            ir = 0
            DO loop = min(kind1,2), 1, -1
              ms = pivot(loop)
              nsc(1) = nsvars
              nr(1) = ir
! Remove chosen variable from linked list (unless it has been
! split).
              IF (ms/=pivt) THEN
                ns = next(ms)
                ls = last(ms)
                next(ms) = 0
                IF (ns>0) last(ns) = ls
                IF (ls>0) THEN
                  next(ls) = ns
                ELSE
                  ipr(count(ms)) = ns
                END IF
              END IF

! Run through the list of the pivotal supervariable, setting
! tree
! pointers and constructing new list of supervariables.
! KP is a pointer to the current position in the old list.
              kp = ipe(ms)
              kp1 = kp + 1
              kp2 = kp + iw(kp)
              DO kp = kp1, kp2
                part = (iw(kp)-1)/n
                ke = iw(kp) - part*n
                IF (flag(ke)/=-1) THEN
                  IF (flag(ke)<=-2) THEN
! KE is an element.
! Link KE or its ancestor to tree if not already linked.
                    ie = ke
                    DO list = 1, n
                      IF (next(ie)==0) THEN
                        GO TO 120
                      ELSE
                        iel = ie
                        ie = -last(ie)
                        IF (ie==me) THEN
                          GO TO 130
                        ELSE
                          last(iel) = -me
                        END IF
                      END IF
                    END DO
120                 next(ie) = -me
                    last(ie) = -me
! Find the relevant part of the list of element KE.
130                 jp = ipe(ke)
                    jp1 = jp + 3
                    jp2 = jp + iw(jp)
                    IF (part/=0) THEN
! Supervariable in full part of element
                      IF (part==2) THEN
! Supervariable in leading part of element
                        jp1 = jp1 + iw(jp+1)
                      ELSE
! Supervariable in trailing part of element
                        jp2 = jp2 - iw(jp+2)
                      END IF
                    END IF
                  ELSE
! We have reached the list of variables
                    jp1 = kp
                    jp2 = kp2
                  END IF

! Search for different supervariables and add them to the
! new
! list.
! This loop is executed once for each element in the list
! and
! once
! for all the supervariables in the list.
                  DO jp = jp1, jp2
                    is = iw(jp)
                    IF (flag(is)>loop) THEN
! IS is not a supervariable or has already been counted.
                      IF (flag(is)==2) THEN
! Supervariable in both rows of the block pivot
                        nvc = nvc + abs(nv(is))
                        nsc(2) = nsc(2) + 1
                        flag(is) = 0
                        count(is) = count(is) - nvpiv
                      ELSE
! New supervariable
                        ir = ir + abs(nv(is))
                        flag(is) = loop
! Store IS in new list.
                        nsvars = nsvars + 1
                        svars(nsvars) = is
! Remove IS from row count linked list
                        ls = last(is)
                        last(is) = 0
                        ns = next(is)
                        next(is) = 0
                        IF (ns>0) last(ns) = ls
                        IF (ls>0) THEN
                          next(ls) = ns
                        ELSE
                          ipr(count(is)) = ns
                        END IF
                        count(is) = count(is) - nvpiv
                      END IF
                    END IF
                  END DO

                  IF (jp2/=kp2 .AND. loop==1) THEN
! See if the pivot is covered by element KE. If so,
! element
! KE will
! be covered by the new element and is no longer needed.
                    IF (kind1/=1) THEN
                      DO jp = jp1, jp2
                        IF (iw(jp)==pivot(2)) GO TO 140
                      END DO
                      GO TO 150
140                   flag(ke) = -1
                      GO TO 160
                    ELSE IF (part==0) THEN
                      flag(ke) = -1
                    END IF
150                 CONTINUE
                  END IF
                END IF
160             CONTINUE
              END DO
            END DO

            IF (kind1==1) THEN
! Complete calculation of counts for a full pivot.
              nsc(2) = nsvars
              nr(1) = ir
              nsc(3) = 0
              nr(3) = 0
            ELSE
! Link PIV(2) into tree.
              next(piv(2)) = -me
              last(piv(2)) = -me
              count(piv(2)) = -count(piv(2))
              IF (kind1==2) THEN
! Complete calculation of counts for a tile pivot.
                nsc(3) = nsvars - nsc(1)
                nr(3) = ir - nr(1)
                nsc(2) = nsc(1)
                nr(2) = nr(1)
                nsc(1) = 0
                nr(1) = ir
              ELSE
! Complete calculation of counts for an oxo pivot.
                nsc(3) = nsvars - nsc(1)
                nr(3) = ir - nr(1) + nvc
                nsc(1) = nsc(1) - nsc(2)
                nr(2) = nr(1)
                nr(1) = ir
! Bring the variables of the leading part to the front.
                k1 = 1
                DO k = 1, nsc(1) + nsc(2)
                  IF (flag(svars(k))==2) THEN
                    ks = svars(k)
                    svars(k) = svars(k1)
                    svars(k1) = ks
                    k1 = k1 + 1
                  END IF
                END DO
              END IF
            END IF

! Deal with the case where the list of the new element has zero
! length.
            IF (nsvars==0) THEN
              ipe(me) = 0
            ELSE
! Run through new list of supervariables looking for elements
! that
! may be absorbed into the new element
! Reduce NFLG by one to cater for these tests.
              IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
              nflg = nflg - 1
              DO k = 1, nsvars
                is = svars(k)
! Run through the list associated with supervariable IS
                kp = ipe(is)
                kp1 = kp + 1
                kp2 = iw(kp) + kp
                DO kp = kp1, kp2
                  part = (iw(kp)-1)/n
                  ke = iw(kp) - n*part
! Exit if we have reached the list of supervariables
                  IF (flag(ke)>=0) THEN
                    GO TO 190
! Cycle if this is a dummy or has already been checked
                  ELSE IF (flag(ke)/=-1) THEN
                    IF (flag(ke)/=-nflg) THEN
                      flag(ke) = -nflg
                      jp = ipe(ke)
                      jp1 = jp + 3
                      jp2 = jp1 + iw(jp+1)
                      jp4 = jp + iw(jp)
                      jp3 = jp4 - iw(jp+2)
                      IF (kind1==1) THEN
! Check for inclusion in a new full element
                        DO jp = jp1, jp4
                          IF (flag(iw(jp))>2) GO TO 180
                        END DO

                      ELSE IF (kind1==2) THEN
! Check for inclusion in a new tile element.
! Check the full part
                        DO jp = jp2, jp3
                          IF (flag(iw(jp))>2) THEN
                            GO TO 180
                          ELSE IF (flag(iw(jp))==1) THEN
                            GO TO 180
                          END IF
                        END DO
                        flg1 = 0
! Check the trailing part
                        DO jp = jp3 + 1, jp4
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            flg1 = flg
                          END IF
                        END DO
! Check the leading part
                        DO jp = jp1, jp2 - 1
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1/=0) GO TO 180
                          END IF
                        END DO
                      ELSE

! Check for inclusion in a new oxo element.
! Check the full part
                        DO jp = jp2, jp3
                          IF (flag(iw(jp))>0) GO TO 180
                        END DO
                        flg1 = 0
! Check the trailing part
                        DO jp = jp3 + 1, jp4
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1==0) flg1 = flg
                            IF (flg/=flg1) GO TO 180
                          END IF
                        END DO
                        flg1 = 3 - flg1
! Check the leading part
                        DO jp = jp1, jp2 - 1
                          flg = flag(iw(jp))
                          IF (flg>2) THEN
                            GO TO 180
                          ELSE IF (flg>0) THEN
                            IF (flg1==3) flg1 = flg
                            IF (flg/=flg1) GO TO 180
                          END IF
                        END DO
                      END IF
! Element KE is absorbed into new element ME.
                      flag(ke) = -1
! Link KE or its ancestor to tree if not already linked.
                      ie = ke
                      DO list = 1, n
                        IF (next(ie)==0) THEN
                          GO TO 170
                        ELSE
                          iel = ie
                          ie = -last(ie)
                          IF (ie==me) THEN
                            GO TO 180
                          ELSE
                            last(iel) = -me
                          END IF
                        END IF
                      END DO
170                   next(ie) = -me
                      last(ie) = -me
                    END IF
                  END IF
180               CONTINUE
                END DO
190             CONTINUE
              END DO

! Run through new list of supervariables revising each
! associated
! list,
! recalculating row counts and removing duplicates.
              DO loop = 1, kind1
! Find the range to be treated and set flags appropriately.
                IF (loop==1) THEN
! Treat middle block
                  k1 = 1 + nsc(1)
                  k2 = k1 + nsc(2) - 1

                ELSE IF (loop==2) THEN
! Treat last block
                  k1 = k2 + 1
                  k2 = nsvars
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = nflg
                  END DO
                ELSE

! Treat first block
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = 1
                  END DO
                  k1 = 1
                  k2 = nsc(1)
                  DO k = k1, k2
                    is = svars(k)
                    IF (nv(is)/=0) flag(is) = nflg
                  END DO
                END IF
! Run through the list of supervariables.
                DO k = k1, k2
                  is = svars(k)
                  bound = count(is) > thresh
                  IF (loop==1) nv(is) = abs(nv(is))
! Reduce NFLG by one to cater for this supervariable.
                  IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
                  nflg = nflg - 1
! Begin with the row count of the new element. Its variables
! must always be counted during the row count calculation
! and they are already flagged with the value 0, 1, or 2.
                  ir = nr(loop)
! Run through the list associated with supervariable IS
                  kp = ipe(is)
! NP points to the next entry in the revised list.
                  np = kp + 1
                  kp1 = kp + 1
                  kp2 = iw(kp) + kp
                  DO kp = kp1, kp2
                    part = (iw(kp)-1)/n
                    ke = iw(kp) - n*part
                    IF (flag(ke)<=-2) THEN
! KE is an element.  Flag it.
                      flag(ke) = -nflg
                      IF ( .NOT. bound) THEN
! Store the part involved in preparation for the test
! for identical
! variables.
                        leaf(ke) = part
! Find the bounds for the search of the element list
                        jp = ipe(ke)
                        jp1 = jp + 3
                        jp2 = jp + iw(jp)
                        IF (part/=0) THEN
! Supervariable in full part of element
                          IF (part==2) THEN
! Supervariable in leading part of element
                            jp1 = jp1 + iw(jp+1)
                          ELSE

! Supervariable in trailing part of element
                            jp2 = jp2 - iw(jp+2)
                          END IF
                        END IF
! Search list of element KE, revising the row count
! when
! new
! variables found.
                        irl = ir
                        DO jp = jp1, jp2
                          js = iw(jp)
! Cycle if JS has been eliminated or already
! counted.
                          IF (flag(js)>nflg) THEN
                            ir = ir + abs(nv(js))
                            flag(js) = nflg
                          END IF
                        END DO
! KP need not be stored in the new list if it is
! linked
! to ME and KE
! made no contribution to the count.
                        IF (ir==irl .AND. last(ke)==-me) GO TO 200
                      END IF
! Store KP in the new list
                      iw(np) = iw(kp)
                      np = np + 1

                    ELSE IF (flag(ke)>=0) THEN
                      GO TO 210
                    END IF
200                 CONTINUE
                  END DO

                  np0 = np
                  GO TO 220
! KE is a supervariable
210               CONTINUE
! Treat the rest of the list associated with supervariable
! IS. It consists entirely of supervariables.
                  np0 = np
                  kp1 = kp
                  DO kp = kp1, kp2
                    ks = iw(kp)
                    IF (flag(ks)>nflg) THEN
! Add to row count, flag supervariable KS and add it to
! new list.
                      ir = ir + abs(nv(ks))
                      flag(ks) = nflg
                      iw(np) = ks
                      np = np + 1
                    END IF
                  END DO
220               IF (bound) ir = count(is)
                  IF (np>kp2) THEN
! List is longer than it was. Copy it to free space.
                    kp = ipe(is)
                    IF (np+iwfr-kp>=lw) THEN
! Compress IW
                      CALL mc68_ma47_compress(n,ipe,flag,iw,iwfr-1,iwfr, &
                        info%n_compressions)
                      IF (np+iwfr-kp>=lw) THEN
! Expand the length of iw
                        szw = lw
                        lw = np + iwfr - kp + 1
                        CALL zb01_resize1(iw,szw,lw,info_zb01)
                        IF (info_zb01%flag<0) THEN
                          IF (info_zb01%flag==-11) THEN
                            info%flag = mc68_err_memory_alloc
                            info%stat = info_zb01%stat
                          ELSE IF (info_zb01%flag==-12) THEN
                            info%flag = mc68_err_memory_dealloc
                            info%stat = info_zb01%stat
                          ELSE IF (info_zb01%flag<-12) THEN
                            info%flag = mc68_err_zb01
                            info%iostat = info_zb01%iostat
                            info%zb01_info = info_zb01%flag
                          END IF
                          RETURN
                        END IF
                      END IF

                      shift = ipe(is) - kp
                      kp = kp + shift
                      kp2 = kp2 + shift
                      np = np + shift
                      np0 = np0 + shift
                    END IF

                    np = np + iwfr - kp
                    np0 = np0 + iwfr - kp
                    ipe(is) = iwfr
                    kp1 = kp
                    DO kp = kp1, kp2
                      iw(iwfr) = iw(kp)
                      iwfr = iwfr + 1
                    END DO
                    iw(iwfr) = 0
                    iwfr = iwfr + 1
                  END IF
! Move first supervariable to end of list, move first
! element
! to end
! of element part of list and add new element to front of
! list.
                  iw(np) = iw(np0)
                  kp = ipe(is)
                  iw(np0) = iw(kp+1)
                  iw(kp+1) = me + (loop-1)*n
! Store the new length of the list.
                  iw(kp) = np - kp
                  IF (ir==0) THEN
! Treat zero row as if it had a nonzero diagonal block
                    ir = -nv(is)
                    nv(is) = ir
                    rank = rank - ir
                    ir = 1
                  END IF

! Unless the Markowitz cost is zero (which means that rows
! not
! previously identical cannot be identical now) or the row
! count is
! above the threshold, check whether row IS is identical to
! another by
! looking in the linked list of supervariables with row
! count
! IR at
! those whose lists have first entry ME and occur in the
! same
! part
! of it. Note that those containing ME come first so the
! search can
! be terminated when a list not starting with ME is found.
                  IF (cmin/=zero) THEN
                    IF (ir<=thresh) THEN
                      js = ipr(ir)
                      DO list = 1, n
                        IF (js<=0) THEN
                          GO TO 260
                        ELSE
                          kp = ipe(js)
                          IF (iw(kp+1)/=me+(loop-1)*n) THEN
                            GO TO 260
                          ELSE
! JS has same row count and is active. Check if
! identical to IS.
                            IF (sign(1,nv(js))==sign(1,nv(is))) THEN
                              kp1 = kp
                              DO kp = kp1 + 2, kp1 + iw(kp1)
                                part = (iw(kp)-1)/n
                                ie = iw(kp) - part*n
! Jump if IE (which may be a
! supervariable or an element)
! is not in the list of IS.
                                IF (abs(flag(ie))>nflg) THEN
                                  GO TO 230
! Jump if JS belongs to another part
! of IE
                                ELSE IF (flag(ie)==-nflg) THEN
                                  IF (part/=leaf(ie)) GO TO 230
                                END IF
                              END DO
                              GO TO 240
                            END IF

230                         js = next(js)
                          END IF
                        END IF
                      END DO
                      GO TO 250

240                   CONTINUE
! Supervariable amalgamation. Row IS is identical to row
! JS.  Regard all variables in the two supervariables as
! being in IS. Set tree pointer, FLAG and NV entries.
250                   ipe(js) = 0
                      nv(is) = nv(is) + nv(js)
                      nv(js) = 0
                      flag(js) = -1
! Replace JS by IS in linked list.
                      ns = next(js)
                      ls = last(js)
                      IF (ns>0) last(ns) = is
                      IF (ls>0) next(ls) = is
                      last(is) = ls
                      next(is) = ns
                      IF (ipr(ir)==js) ipr(ir) = is
                      count(is) = ir
                      next(js) = -leaf(is)
                      leaf(is) = leaf(js)
                      last(js) = -is
                      GO TO 270
                    END IF
                  END IF
! Insert IS into linked list of supervariables of same row
! count.
260               ns = ipr(ir)
                  IF (ns>0) last(ns) = is
                  next(is) = ns
                  ipr(ir) = is
                  last(is) = 0
                  minr = min(minr,ir)
                  IF (nv(is)>0) minf = min(minf,ir)
                  count(is) = ir
270               CONTINUE
                END DO
              END DO

! Reset flags for supervariables in newly created element,
! remove those absorbed into others, and store the new list.
              IF (iwfr+nsvars+3>=lw) THEN
! Compress IW
                CALL mc68_ma47_compress(n,ipe,flag,iw,iwfr-1,iwfr, &
                  info%n_compressions)
                IF (iwfr+nsvars+3>=lw) THEN
                  szw = lw
                  lw = iwfr + nsvars + 4
                  CALL zb01_resize1(iw,szw,lw,info_zb01)
                  IF (info_zb01%flag<0) THEN
                    IF (info_zb01%flag==-11) THEN
                      info%flag = mc68_err_memory_alloc
                      info%stat = info_zb01%stat
                    ELSE IF (info_zb01%flag==-12) THEN
                      info%flag = mc68_err_memory_dealloc
                      info%stat = info_zb01%stat
                    ELSE IF (info_zb01%flag<-12) THEN
                      info%flag = mc68_err_zb01
                      info%iostat = info_zb01%iostat
                      info%zb01_info = info_zb01%flag
                    END IF
                    RETURN
                  END IF
! GO TO 280
                ELSE
                END IF
              END IF

              ip = iwfr
              iwfr = iwfr + 3
              k2 = 0
              DO loop = 1, 3
                k1 = k2 + 1
                k2 = k1 + nsc(loop) - 1
                DO k = k1, k2
                  is = svars(k)
                  IF (nv(is)==0) THEN
                    nsc(loop) = nsc(loop) - 1
                  ELSE
                    flag(is) = nflg
                    iw(iwfr) = is
                    iwfr = iwfr + 1
                  END IF
                END DO
              END DO
              iw(ip) = iwfr - ip - 1
              iw(ip+1) = nsc(1)
              iw(ip+2) = nsc(3)
              flag(me) = -nflg
! Set pointer for new element.
              ipe(me) = ip
            END IF
          END IF
        END DO

280     IF (rank<n) info%flag = info%flag + mc68_warn_rank
        info%n_zero_eigs = n - rank


        DEALLOCATE (flag,ipr,last,leaf,svars,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

      END SUBROUTINE mc68_ma47_analyse

      SUBROUTINE mc68_ma47_treesearch(n,father,count,ne,perm,nemin,info)
! Given the tree provided by MC68_MA47_ANALYSE/KD, which has a node
! for every
! variable, perform depth-first search to find pivot order and
! construct a new tree whose nodes correspond to block eliminations,
! and whose depth-first-search post ordering is 1,2,...,NODES.

! N must be set to the matrix order. It is not altered.
! FATHER(I) must be set equal to -(father of node I) in the given tree
! or zero if node I is a root. It is altered to point to its next
! younger brother if it has one. Finally, it is set to the father
! of node I in the new tree, or NODES+1 if node I is a root.
! COUNT(I) must be set to the row count at the time of elimination if
! I
! is a node at which an elimination takes place, negated for the
! second node of a tile or oxo pivot. It is revised for a node at
! which an amalgamaton occurs, but is otherwise unchanged.
! NE(I) must be set to the number of variables eliminated at node I,
! negated for a defective variable. It is revised for a node at
! which an amalgamaton occurs, but is otherwise unchanged.
! PERM need not be set.

! Local variables
! SON(I) need not be set. It is used to hold
! (eldest son of node I) if it has one and 0 otherwise.
! NODE(I) need not be set. It is used temporarily to hold the row
! count of the secondary variable of a tile or oxo , -1 for a
! node merged into its father and 0 otherwise. It is finally
! set to hold the elimination stage of variable I, negated for a
! defective variable.
! NA need not be set. NA(I) is used temporarily to hold the next
! younger
! brother of node I if it has one or -(father) if it does not.
! NA(STAGE) and NA(LEVEL) are used to hold the numbers of elements
! assembled at stage STAGE and level LEVEL of the elimination.
! MARK need not be set. The indices of the root nodes are stored in
! MARK(I), I=NR,N. MARK(STAGE) eventually is set to the Markowitz
! cost at stage STAGE of the elimination.
! NODES need not be set. It is set to the number of elimination steps.
! ICNTL must be set by the user as follows and is not altered. The
! only entry used is
! ICNTL(6) Two nodes of the assembly tree are merged only if both
! involve less than ICNTL(6) eliminations.
! COUNT2 Row count of the second row of a tile or oxo.
! I   Index of current tree node.
! IBRTHR  Brother of node I.
! IFATHR  Father of node I.
! ISON Son of node I.
! J Temporary variable.
! K Position in elimination sequence
! L Temporary variable.
! LDIAG Control for amount of information output.
! LEVEL holds the current tree level. Roots are at level N, their sons
! are at level N-1, etc.
! MP Stream number for warnings and diagnostics.
! NDE Do index when looping through new tree.
! NEMIN is used to control the amalgamation process between
! a son and its father (if the number of fully summed
! variables of both nodes is smaller than NEMIN).
! NST Number of structured pivots in the active chain of the
! depth-first search.
! NR Position of next root.
! NRL The indices of the root nodes are stored in MARK(I), I=NRL,N.
! NVPIV holds the number of variables in the pivot (each half in the
! case of an oxo).
! PERM1(LEVEL) is used temporarily to hold the node at level LEVEL of
! the
! active chain. PERM1(K) is eventually set to the variable that is in
! position K of the elimination.
! STAGE holds the current elimination stage.
! Set local print variables
! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n, nemin
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: count(n), father(n), ne(n)
        INTEGER, INTENT (OUT) :: perm(n)
        TYPE (mc68_info), INTENT (INOUT) :: info
! ..
! .. Local Scalars ..
        INTEGER :: count2, i, ibrthr, ifathr, ison, j, k, l, ldiag, level, mp, &
          nr, nrl, nst, nvpiv, stage
        INTEGER, ALLOCATABLE, DIMENSION (:) :: mark, na, node, son, perm1
! ..
        mp = -1
        ldiag = 0
        IF (mp<=0) ldiag = 0


        ALLOCATE (mark(n),na(n),node(n),son(n),perm1(n),STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF


! Set SON(I) to (eldest son of node I) and NA(I) to next younger
! brother of node I if it has one or -(father) if it does not.
! Find the root nodes.
! Load all the eliminations of any oxo and tile pivot on its first
! node and store the second count in NODE.
        DO i = 1, n
          son(i) = 0
          node(i) = 0
        END DO
        nrl = n + 1
        DO i = 1, n
          ifathr = -father(i)
          na(i) = -ifathr
          IF (ifathr==0) THEN
! We have a root
            nrl = nrl - 1
            mark(nrl) = i
          ELSE
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) na(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
            IF (count(i)<0) THEN
! Second node of tile or oxo. Load the eliminations on the first
! node
! and store the second count.
              ne(ifathr) = ne(ifathr)*2
              ne(i) = 0
              node(ifathr) = -count(i)
            END IF
          END IF
        END DO

! Depth-first search looking for father-son pairs that can be merged.
! Adjust NE, COUNT, and NODE to correspond.
! PERM1(LEVEL) is used temporarily to hold the node at level
! LEVEL of the active chain.
        nr = nrl
        i = 0
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            nst = 0
            perm1(level) = i
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used.
          DO l = 1, n
            IF (son(i)<=0) THEN
              GO TO 10
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              IF (node(i)/=0) nst = nst + 1
              level = level - 1
              perm1(level) = i
            END IF
          END DO
! Jump if node has no eliminations
10        nvpiv = ne(i)
          IF (nvpiv/=0) THEN
! Jump if at a root.
            IF (level/=n) THEN
              ifathr = perm1(level+1)
              IF (node(i)/=0 .OR. node(ifathr)/=0) THEN
                IF (node(i)==0 .OR. node(ifathr)==0) THEN
                  GO TO 20
                ELSE
                  IF (nvpiv>0 .AND. ne(ifathr)>0) THEN
! Both nodes are tiles.
! Do not merge if this would cause fill.
                    IF (node(i)-nvpiv/2/=node(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
                      GO TO 20
                    END IF
                  ELSE IF (nvpiv<0 .AND. ne(ifathr)<0) THEN
! Both nodes are oxos.
                    nvpiv = -nvpiv/2
! Do not merge unless sure that this would not cause fill.
                    IF (node(i)-nvpiv/=node(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
                      GO TO 20
                    ELSE IF (count(i)==node(i)) THEN
                      GO TO 20
                    END IF
                  ELSE
                    GO TO 20
                  END IF
! Merge two tiles or two oxos.
                  node(ifathr) = node(i)
                  nst = nst - 1
                END IF
! Both nodes are full
! Merge the nodes if this would cause no more fill
              ELSE IF (count(i)-nvpiv/=count(ifathr)) THEN
! Do not merge the nodes if there is a structured pivot on
! path
! to root
                IF (nst>0) THEN
                  GO TO 20
! Do not merge the nodes if either node has more than NEMIN
! eliminations
                ELSE IF (nvpiv>=nemin) THEN
                  GO TO 20
                ELSE IF (ne(ifathr)>=nemin) THEN
                  GO TO 20
                END IF
              END IF
! Merge the two nodes
              ne(ifathr) = ne(ifathr) + ne(i)
              ne(i) = 0
              IF (ldiag>4) WRITE (mp,fmt='(A,2I5)') ' Merging nodes', i, &
                ifathr
              count(ifathr) = count(ifathr) + nvpiv
              node(i) = -1
            END IF
          END IF
! Go to next brother, or failing this to the father.
20        ibrthr = na(i)
          IF (node(i)>0) nst = nst - 1
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            perm1(level) = ibrthr
            i = ibrthr
            IF (node(i)>0) nst = nst + 1
          ELSE
! Go to father of node I
            level = level + 1
            i = -ibrthr
          END IF
        END DO
        DO i = 1, n
          son(i) = 0
        END DO

! Set SON(I) to (eldest son of node I) and NA(I) to next younger
! brother of node I if it has one or -(father) if it does not.
        DO i = 1, n
          ifathr = -father(i)
          na(i) = -ifathr
          IF (ifathr/=0) THEN
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) na(i) = ibrthr
! Set father-to-son pointer
            son(ifathr) = i
          END IF
        END DO

! Depth-first search in which FATHER is revised to accord with the
! merges. PERM1(LEVEL) is used temporarily to hold the node at level
! LEVEL of the active chain.
        i = 0
        nr = nrl
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            perm1(n) = i
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used.
          DO l = 1, n
            IF (son(i)<=0) THEN
              GO TO 30
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              father(i) = -perm1(level)
              IF (node(i)>=0) THEN
                level = level - 1
                perm1(level) = i
              END IF
            END IF
          END DO

! Go to next brother, or failing this to the father.
30        ibrthr = na(i)
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            IF (node(i)<0) level = level - 1
            i = ibrthr
            perm1(level) = i
            father(i) = -perm1(level+1)
            IF (node(i)<0) level = level + 1
          ELSE
! Go to father of node I
            IF (node(i)>=0) level = level + 1
            i = -ibrthr
          END IF
        END DO

! Set SON(I) to (eldest son of node I) and FATHER(I) to next younger
! brother of node I if it has one.
        DO i = 1, n
          son(i) = 0
        END DO
! First pass is for nodes without eliminations.
        DO i = 1, n
          IF (ne(i)==0 .AND. count(i)>=0) THEN
            ifathr = -father(i)
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) father(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
          END IF
        END DO

! Second pass is for second nodes of tile and oxo pivots.
        DO i = 1, n
          IF (count(i)<0) THEN
            ifathr = -father(i)
            ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
            IF (ibrthr>0) father(i) = ibrthr
! Set father-to-eldest-son pointer
            son(ifathr) = i
          END IF
        END DO

! Third pass is for elimination nodes.
        DO i = 1, n
          IF (ne(i)/=0) THEN
            ifathr = -father(i)
            IF (ifathr/=0) THEN
              ibrthr = son(ifathr)
! If node I has a next younger brother, set pointer to him.
              IF (ibrthr>0) father(i) = ibrthr
! Set father-to-son pointer
              son(ifathr) = i
            END IF
          END IF
        END DO

! Depth-first search.
! STAGE holds the current elimination stage. The number
! of assemblies for nodes of the active chain is accumulated
! temporarily in NA(LEVEL), for tree level N, N-1,..., and is
! transfered to NA(STAGE) when we reach the appropriate stage STAGE.
        stage = 1
! I is the current node.
        i = 0
        nr = nrl
        DO k = 1, n
          IF (i<=0) THEN
! Pick up next root.
            i = mark(nr)
            nr = nr + 1
            level = n
            na(n) = 0
          END IF
! Go to son for as long as possible, clearing father-son pointers
! in SON as each is used and setting NA(LEVEL)=0 for all levels
! reached.
          l = level
          DO level = l, 1, -1
            IF (son(i)<=0) THEN
              GO TO 40
            ELSE
              ison = son(i)
              son(i) = 0
              i = ison
              na(level-1) = 0
            END IF
          END DO
! Record variable in position I in the order.
40        perm1(k) = i
          count2 = node(i)
          node(i) = stage
! Jump if node has no eliminations
          IF (ne(i)/=0) THEN
            IF (level<n) na(level+1) = na(level+1) + 1
            na(stage) = na(level)
            IF (count2==0) THEN
! Full pivot
              mark(stage) = (count(i)-1)**2
            ELSE IF (ne(i)>0) THEN
! Tile pivot
              DO j = k - ne(i) + 1, k - ne(i)/2
                node(perm1(j)) = -node(perm1(j))
              END DO
              mark(stage) = (count(i)+count2-3)*(count2-1)
            ELSE
! Oxo pivot
              DO j = k + ne(i) + 1, k
                node(perm1(j)) = -node(perm1(j))
              END DO
              mark(stage) = (count(i)-1)*(count2-1)
            END IF
            stage = stage + 1
          END IF
! Go to next brother, or failing this to the father.
          ibrthr = father(i)
          IF (ibrthr>0) THEN
! Node I has a younger brother. Go to him.
            na(level) = 0
            i = ibrthr
          ELSE
! Go to father of node I
            level = level + 1
            ison = i
            i = -ibrthr
          END IF
        END DO

! Search for 2x2 pivots
        i = 1
        DO WHILE (i<n)
          IF (node(perm1(i))>0) THEN
            i = i + 1
          ELSE
            perm1(i) = -perm1(i)
            perm1(i+1) = -perm1(i+1)
            i = i + 2
          END IF
        END DO
        DO i = 1, n
          IF (perm1(i)>0) THEN
            perm(perm1(i)) = i
          ELSE
            perm(-perm1(i)) = -i
          END IF
        END DO

        DEALLOCATE (mark,na,node,son,perm1,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF


      END SUBROUTINE mc68_ma47_treesearch

      SUBROUTINE mc68_ma47_merge(n,ipe,iw,lw,nv,next,last,leaf,flag,var,svar)

! Merge variables whose rows are identical into supervariables. This
! is done by updating the supervariable structure for the submatrix
! of the first J-1 columns to that for the first J columns,
! J=1,2,...,N. If a nontrivial supervariable IS is involved in
! column J, all of its variables that are involved in column J are
! removed to make a new supervariable.

! N must be set to the matrix order. It is not altered.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N. On return, it is unchanged except that if I is a variable
! that has been absorbed into a supervariable, IPE(I) = 0.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. There must be no duplicate indices
! in a list. IW is not altered.
! LW must be set to the length of IW. It is not altered.
! NV must be set on entry so that NV(I) = 1 for a nondefective
! variable
! (nonzero diagonal entry) and NV(I) = -1 otherwise, I = 1,..., N.
! On return:
! if I is a variable that has been absorbed, NV(I) = 0;
! if I is a supervariable, NV(I) holds its number of variables,
! negated for a defective supervariable.
! NEXT need not be set. During execution, NEXT(I) is the next
! variable in a circular list of variables in a supervariable.
! On return, if variable I has been absorbed and J is the next
! variable in its supervariable, NEXT(I)=-J;
! LAST need not be set. During execution, LAST(I) is the previous
! variable in a circular list of variables in a supervariable.
! On return, if variable I has been absorbed and J is the next
! variable in its supervariable, LAST(I)=-J;
! LEAF need not be set. On return, if IS is a supervariable, LEAF(IS)
! holds its first variable.
! FLAG is used as workspace for supervariable flags. If IS has already
! been encountered in column J, FLAG(IS) = J.
! VAR is used as workspace. If IS is a supervariable involved in
! column
! J, VAR(IS) is the first variable to be removed from supervariable
! IS. If IS is not a supervariable, VAR(IS) is the next free
! supervariable index in a chain of such indices.
! SVAR is used as workspace. SVAR(I) is the supervariable to which I
! belongs.

! Local variables

! FREE The first free supervariable index (head of chain)
! I    Row index.
! IS   Supervariable.
! JS   Supervariable.
! J    Column index.
! K    DO index.
! KK   Temporary variable.
! LS   Last (previous) variable in circular list.
! NS   Next variable in circular list.

! Begin by setting SVAR, LAST, and NEXT to represent all variables
! belonging to supervariable 1. Also initialize FLAG and set
! FREE and VAR as a chain of indices that are free for use as
! supervariable names.
! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: iw(lw)
        INTEGER, INTENT (INOUT) :: ipe(n), nv(n)
        INTEGER, INTENT (OUT) :: flag(n), last(n), leaf(n), next(n), svar(n), &
          var(n)
! ..
! .. Local Scalars ..
        INTEGER :: free, i, is, j, js, k, kk, ls, ns
! ..
! .. Intrinsic Functions ..
        INTRINSIC sign
! ..
        DO i = 1, n
          svar(i) = 1
          last(i) = i - 1
          next(i) = i + 1
          flag(i) = 0
          var(i) = i + 1
        END DO
        last(1) = n
        next(n) = 1
        free = 2

! Scan the columns in turn, splitting the supervariables that are
! involved.
        DO j = 1, n
          kk = ipe(j)
          DO k = kk + 1, kk + iw(kk)
            i = iw(k)
            is = svar(i)
            IF (flag(is)/=j) THEN
! First occurrence of supervariable IS for column J
              flag(is) = j
              IF (next(i)/=i) THEN
! No action needed since IS has I as its only variable.
! Establish new supervariable
                js = free
                free = var(js)
                var(is) = i
                svar(i) = js
! Remove I from old circular list
                ns = next(i)
                ls = last(i)
                next(ls) = ns
                last(ns) = ls
! Make new circular list
                next(i) = i
                last(i) = i
              END IF
            ELSE
! Subsequent occurrence of IS for column J
! Remove I from old circular list
              IF (next(i)==i) THEN
! Supervariable now empty
                ns = var(is)
                var(is) = free
                free = is
              ELSE
                ns = next(i)
                ls = last(i)
                next(ls) = ns
                last(ns) = ls
                ns = var(is)
              END IF
! Add I to new list
              ls = last(ns)
              next(ls) = i
              next(i) = ns
              last(ns) = i
              last(i) = ls
              svar(i) = svar(ns)
            END IF
          END DO
        END DO

! Set the data in final format.
        DO is = 1, n
          leaf(is) = is
! No action for a trivial supervariable
          IF (last(is)/=is) THEN
! No action if already treated
            IF (last(is)>=0) THEN
! IS is a nontrivial supervariable
              ls = last(is)
              leaf(is) = ls
              DO k = 1, n
                i = ls
                IF (i==is) THEN
                  GO TO 10
                ELSE
                  ipe(i) = 0
                  ls = last(i)
                  next(i) = -ls
                  last(i) = -ls
                  nv(i) = 0
                END IF
              END DO
10            nv(is) = sign(k,nv(is))
            END IF
          END IF
        END DO
      END SUBROUTINE mc68_ma47_merge

      SUBROUTINE mc68_ma47_reset(n,flag,nflg)
! Reset flag values to +/- N*3.
! N Matrix order. Unchanged.
! FLAG Unchaged if FLAG(I) = -1, 0, 1, or 2.
! Otherwise, reset to +/- N*3.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: n
        INTEGER, INTENT (OUT) :: nflg
! ..
! .. Array Arguments ..
        INTEGER, INTENT (INOUT) :: flag(n)
! ..
! .. Local Scalars ..
        INTEGER :: i, n3
! ..
        n3 = n*3
        DO i = 1, n
          IF (flag(i)>2) flag(i) = n3
          IF (flag(i)<=-2) flag(i) = -n3
        END DO
        nflg = n3
      END SUBROUTINE mc68_ma47_reset

      SUBROUTINE mc68_ma47_thresh(thresh,newthr,n,ipe,iw,lw,count,nv,next, &
          last,ipr,flag,nflg)
! Alter the threshold for calculation of row counts.

! THRESH Old threshold. Changed to NEWTHR.
! NEWTHR New threshold. It is not altered.
! N must be set to the matrix order. It is not altered.
! IPE(I) must be set to the position in IW of the list for row I,
! I=1,N.
! It is not altered.
! IW must be set on entry to hold lists of entries by rows, each list
! being headed by its length. It is not altered.
! LW must be set to the length of IW. It is not altered.
! COUNT is used to hold the row counts in the reduced matrix.
! COUNT(IS) is recalculated if THRESH < COUNT(IS) <= NEWTHR.
! NV: if JS is a supervariable, NV(JS) holds its number of variables,
! negated for a defective supervariable. It is not altered.
! NEXT: if JS is a supervariable that has not been eliminated or
! absorbed, NEXT(JS) is the next such supervariable having the
! same row count, or zero if JS is last in its list. Revised
! if a recalculated row count makes it necessary.
! LAST: if JS is a supervariable that has not been eliminated or
! absorbed, LAST(JS) is the previous such supervariable having
! the same row count, or zero if JS is last in its list.
! Revised if a recalculated row count makes it necessary.
! IPR: IPR(IR) is the first supervariable with row count IR or zero
! if there are none. Revised if a recalculated row count makes it
! necessary.
! FLAG:  if IS is a supervariable, FLAG(IS) >= 0;
! if IE is an element, FLAG(IE) < -1;
! if I is neither a supervariable nor an element, FLAG(I) = -1.
! The values of the supervariable flags may change while remaining
! positive. Otherwise, unchanged.
! NFLG is used for the current flag value in array FLAG, as in
! MC68_MA47_ANALYSE.
! Local variables

! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: lw
        INTEGER, INTENT (IN) :: n, newthr
        INTEGER, INTENT (INOUT) :: nflg, thresh
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: ipe(n), iw(lw), nv(n)
        INTEGER, INTENT (INOUT) :: count(n), flag(n), ipr(n), last(n), next(n)
! ..
! .. Local Scalars ..
        INTEGER :: ir, is, jp, jp1, jp2, k, ke, kp, kp2, ls, ms, ns, part
! ..
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs
! ..
        DO ms = 1, n
          IF (flag(ms)>=0) THEN
            IF (count(ms)>thresh) THEN
              IF (count(ms)<=newthr) THEN
! Find the row count afresh
                ir = 0
                IF (nflg<=4) CALL mc68_ma47_reset(n,flag,nflg)
! Reduce NFLG by one to cater for this supervariable.
                nflg = nflg - 1
                k = ipe(ms)
                kp2 = k + iw(k)
                DO kp = k + 1, kp2
                  part = (iw(kp)-1)/n
                  ke = iw(kp) - part*n
                  IF (flag(ke)/=-1) THEN
                    IF (flag(ke)<=-2) THEN
! KE is an element.
                      jp = ipe(ke)
                      jp1 = jp + 3
                      jp2 = jp + iw(jp)
                      IF (part/=0) THEN
! Supervariable in full part of element
                        IF (part==2) THEN
! Supervariable in leading part of element
                          jp1 = jp1 + iw(jp+1)
                        ELSE
! Supervariable in trailing part of element
                          jp2 = jp2 - iw(jp+2)
                        END IF
                      END IF
                    ELSE
! We have reached the list of variables
                      jp1 = kp
                      jp2 = kp2
                    END IF
! Search the variable list.
                    DO jp = jp1, jp2
                      is = iw(jp)
                      IF (flag(is)>nflg) THEN
                        flag(is) = nflg
                        ir = ir + abs(nv(is))
                      END IF
                    END DO
                    IF (jp2==kp2 .OR. ir>newthr) GO TO 10
                  END IF
                END DO
10              IF (ir/=count(ms)) THEN
! Remove MS from linked list
                  ns = next(ms)
                  ls = last(ms)
                  next(ms) = 0
                  IF (ns>0) last(ns) = ls
                  IF (ls>0) THEN
                    next(ls) = ns
                  ELSE
                    ipr(count(ms)) = ns
                  END IF
! Insert MS into linked list of supervariables of same row
! count.
                  ns = ipr(ir)
                  IF (ns>0) last(ns) = ms
                  next(ms) = ns
                  ipr(ir) = ms
                  last(ms) = 0
                  count(ms) = ir
                END IF
              END IF
            END IF
          END IF
        END DO
        thresh = newthr
      END SUBROUTINE mc68_ma47_thresh

      SUBROUTINE mc68_ma47_compress(n,ipe,flag,iw,lw,iwfr,ncmpa)
! Compress lists held by MC68_MA47_ANALYSE and MA47KD in IW, removing
! inactive
! variables from element lists. Adjust pointers in IPE to correspond.
! N is the matrix order. It is not altered.
! IPE(I) points to the position in IW of the start of list I or is
! zero if there is no list I. On exit it points to the new position.
! FLAG holds element and supervariable flags:
! if IS is an active supervariable, FLAG(IS) >= 0;
! if IS is inactive, FLAG(IS) = -1;
! if IE is an element, FLAG(IE) < -1;
! It is not altered.
! IW holds the lists, each headed by its length. Every entry must be
! nonnegative. On output, the same lists are held, but they are
! now compressed together.
! LW holds the length of IW. It is not altered.
! IWFR need not be set on entry. On exit it points to the first free
! location in IW.

! Local variables
! I    Row index.
! IR   DO index.
! K    Temporary variable.
! L    Temporary variable.
! LEN1 Length of the leading zero part of element list
! LEN2 Length of the middle part of element list
! LEN3 Length of the trailing zero part of element list
! LWFR points just beyond the end of the uncompressed file.

! .. Scalar Arguments ..
        INTEGER, INTENT (IN) :: lw, n
        INTEGER, INTENT (INOUT) :: ncmpa
        INTEGER, INTENT (OUT) :: iwfr
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: flag(n)
        INTEGER, INTENT (INOUT) :: ipe(n), iw(lw)
! ..
! .. Local Scalars ..
        INTEGER :: i, ir, k, l, len1, len2, len3, lwfr
! ..
        ncmpa = ncmpa + 1
! Prepare for compressing by storing the lengths of the lists in IPE
! and
! setting the first entry of each list to -(list number).
        DO i = 1, n
          l = ipe(i)
          IF (l>0 .AND. flag(i)/=-1) THEN
            ipe(i) = iw(l)
            iw(l) = -i
          END IF
        END DO

! Compress
! IWFR points just beyond the end of the compressed file.
        iwfr = 1
        lwfr = 1
        DO ir = 1, n
! Search for the next negative entry.
          DO k = lwfr, lw
            IF (iw(k)<0) GO TO 10
          END DO
          RETURN
! Pick up entry number, store length in new position, set new
! pointer
! and prepare to copy list.
10        i = -iw(k)
          iw(iwfr) = ipe(i)
          ipe(i) = iwfr
          iwfr = iwfr + 1
          IF (flag(i)<=-2) THEN
! We have an element list. Remove dummy entries.
            l = iwfr - 1
            len1 = iw(k+1)
            iw(l+1) = len1
            len3 = iw(k+2)
            iw(l+2) = len3
            len2 = iw(l) - len1 - len3 - 2
            iwfr = l + 3
            DO lwfr = k + 3, k + 2 + len1
              IF (flag(iw(lwfr))<0) THEN
                iw(l+1) = iw(l+1) - 1
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO

            k = lwfr
            DO lwfr = k, k - 1 + len2
              IF (flag(iw(lwfr))<0) THEN
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO

            k = lwfr
            DO lwfr = k, k - 1 + len3
              IF (flag(iw(lwfr))<0) THEN
                iw(l+2) = iw(l+2) - 1
                iw(l) = iw(l) - 1
              ELSE
                iw(iwfr) = iw(lwfr)
                iwfr = iwfr + 1
              END IF
            END DO
          ELSE
! Copy list to new position.
            DO lwfr = k + 1, k + iw(iwfr-1)
              iw(iwfr) = iw(lwfr)
              iwfr = iwfr + 1
            END DO
          END IF
        END DO

      END SUBROUTINE mc68_ma47_compress

      SUBROUTINE amdd(n,iwlen,pe,pfree,len,iw,elen,icntl,info)

! -------------------------------------------------------------------
! AMDD is a modified version of MC47 version 1. Dense rows are detected
! and removed. AMDD is then applied to the remaining matrix. The dense rows are
! then appended to the end of this ordering.

! We use the term Le to denote the set of all supervariables in element
! E.
! A row is declared as dense if removing it from the current matrix
! will result in a significant reduction in the mean degree of the remaining
! rows
! A row is sparse if it is not dense.
! -------------------------------------------------------------------

! .. Parameters ..
        INTEGER, PARAMETER :: wp = kind(0.0D0)
        REAL (wp) :: zero
        PARAMETER (zero=0E0_wp)
! ..
! .. Scalar Arguments ..
        INTEGER (long), INTENT (IN) :: iwlen
        INTEGER, INTENT (IN) :: n
        INTEGER, INTENT (INOUT) :: pfree
! ..
! .. Array Arguments ..
        INTEGER, INTENT (IN) :: icntl(10)
        INTEGER, INTENT (INOUT) :: iw(iwlen), len(n), pe(n)
        INTEGER, INTENT (OUT) :: elen(n)

        TYPE (mc68_info), INTENT (INOUT) :: info

! N must be set to the matrix order.
! Restriction:  N .ge. 1

! IWLEN must be set to the length of IW. It is not altered. On input,
! the matrix is stored in IW (1..PFREE-1).
! *** We do not recommend running this algorithm with ***
! ***      IWLEN .LT. PFREE + N.                      ***
! *** Better performance will be obtained if          ***
! ***      IWLEN .GE. PFREE + N                       ***
! *** or better yet                                   ***
! ***      IWLEN .GT. 1.2 * PFREE                     ***
! Restriction: IWLEN .GE. PFREE-1

! PE(i) must be set to the the index in IW of the start of row I, or be
! zero if row I has no off-diagonal entries. During execution,
! it is used for both supervariables and elements:
! * Principal supervariable I:  index into IW of the
! list of supervariable I.  A supervariable
! represents one or more rows of the matrix
! with identical pattern.
! * Non-principal supervariable I:  if I has been absorbed
! into another supervariable J, then PE(I) = -J.
! That is, J has the same pattern as I.
! Note that J might later be absorbed into another
! supervariable J2, in which case PE(I) is still -J,
! and PE(J) = -J2.
! * Unabsorbed element E:  the index into IW of the list
! of element E.  Element E is created when
! the supervariable of the same name is selected as
! the pivot.
! * Absorbed element E:  if element E is absorbed into element
! E2, then PE(E) = -E2.  This occurs when one of its
! variables is eliminated and when the pattern of
! E (that is, Le) is found to be a subset of the pattern
! of E2 (that is, Le2).  If element E is "null" (it has
! no entries outside its pivot block), then PE(E) = 0.

! On output, PE holds the assembly tree/forest, which implicitly
! represents a pivot order with identical fill-in as the actual
! order (via a depth-first search of the tree). If NV(I) .GT. 0,
! then I represents a node in the assembly tree, and the parent of
! I is -PE(I), or zero if I is a root. If NV(I)=0, then (I,-PE(I))
! represents an edge in a subtree, the root of which is a node in
! the assembly tree.

! PFREE must be set to the position in IW of the first free variable.
! During execution, additional data is placed in IW, and PFREE is
! modified so that components  of IW from PFREE are free.
! On output, PFREE is set equal to the size of IW that would have
! caused no compressions to occur.  If NCMPA is zero, then
! PFREE (on output) is less than or equal to IWLEN, and the space
! IW(PFREE+1 ... IWLEN) was not used. Otherwise, PFREE (on output)
! is greater than IWLEN, and all the memory in IW was used.

! LEN(I) must be set to hold the number of entries in row I of the
! matrix, excluding the diagonal.  The contents of LEN(1..N) are
! undefined on output.

! IW(1..PFREE-1) must be set to  hold the patterns of the rows of
! the matrix.  The matrix must be symmetric, and both upper and
! lower triangular parts must be present.  The diagonal must not be
! present.  Row I is held as follows:
! IW(PE(I)...PE(I) + LEN(I) - 1) must hold the list of
! column indices for entries in row I (simple
! supervariables), excluding the diagonal.  All
! supervariables start with one row/column each
! (supervariable I is just row I). If LEN(I) is zero on
! input, then PE(I) is ignored on input. Note that the
! rows need not be in any particular order, and there may
! be empty space between the rows.
! During execution, the supervariable I experiences fill-in. This
! is represented by constructing a list of the elements that cause
! fill-in in supervariable I:
! IE(PE(i)...PE(I) + ELEN(I) - 1) is the list of elements
! that contain I. This list is kept short by removing
! absorbed elements. IW(PE(I)+ELEN(I)...PE(I)+LEN(I)-1)
! is the list of supervariables in I. This list is kept
! short by removing nonprincipal variables, and any entry
! J that is also contained in at least one of the
! elements in the list for I.
! When supervariable I is selected as pivot, we create an element E
! of the same name (E=I):
! IE(PE(E)..PE(E)+LEN(E)-1) is the list of supervariables
! in element E.
! An element represents the fill-in that occurs when supervariable
! I is selected as pivot.
! CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
! The contents of IW are undefined on output.

! ELEN(I) need not be set. See the description of IW above. At the
! start of execution, ELEN(I) is set to zero. For a supervariable,
! ELEN(I) is the number of elements in the list for supervariable
! I. For an element, ELEN(E) is the negation of the position in the
! pivot sequence of the supervariable that generated it. ELEN(I)=0
! if I is nonprincipal.
! On output ELEN(1..N) holds the inverse permutation (the same
! as the 'INVP' argument in Sparspak). That is, if K = ELEN(I),
! then row I is the Kth pivot row.  Row I of A appears as the
! (ELEN(I))-th row in the permuted matrix, PAP^T.

! ICNTL is an INTEGER array of length 10 that contains control
! parameters and must be set by the user. 
!     ICNTL(1) - ICNTL(3) are not used.

!     ICNTL(4) controls the choice of AMD algorithm
!   <= 0 No checking for dense rows performed
!   > 0 Corresponds to automatic setting of the minimum density
!     requirement. mc68_order sets this to 1.

!     ICNTL(5) defines the largest positive
!     integer that your computer can represent (-iovflo should also
!     be representable). HUGE(1) in Fortran 95. 

! Local arrays:
! ---------------

! NV(I) During execution, ABS(NV(I)) is equal to the
! number of rows represented by the principal supervariable I. If I
! is a nonprincipal variable, then NV(I) = 0. Initially, NV(I) = 1
! for all I.  NV(I) .LT. 0 signifies that I is a principal variable
! in the pattern Lme of the current pivot element ME. On termination,
! NV(E) holds the true degree of element E at the time it was
! created (including the diagonal part).

! LAST(I) In a degree list, LAST(I) is the
! supervariable preceding I, or zero if I is the head of the list.
! In a hash bucket, LAST(I) is the hash key for I. LAST(HEAD(HASH))
! is also used as the head of a hash bucket if HEAD(HASH) contains
! a degree list (see HEAD, below).
! On output, LAST(1..N) holds the permutation (the same as the
! 'PERM' argument in Sparspak). That is, if I = LAST(K), then row I
! is the Kth pivot row.  Row LAST(K) of A is the K-th row in the
! permuted matrix, PAP^T.


! DEGREE If I is a supervariable and sparse,
! then DEGREE(I) holds the current approximation of the external
! degree of row I (an upper bound). The external degree is the
! number of entries in row I, minus ABS(NV(I)) (the diagonal
! part). The bound is equal to the external degree if ELEN(I) is
! less than or equal to two. We also use the term "external degree"
! for elements E to refer to |Le \ Lme|. If I is full in the reduced
! matrix, then DEGREE(I)=N+1. If I is dense in the reduced matrix,
! then DEGREE(I)=N+1+last_approximate_external_deg of I.
! All dense rows are stored in the list pointed by HEAD(N).
! Quasi dense rows are stored first, and are followed by full rows
! in the reduced matrix. LASTD holds the last row in
! this list of dense rows or is zero if the list is empty.

! HEAD(DEG) is used for degree lists.
! HEAD(DEG) is the first supervariable in a degree list (all
! supervariables I in a degree list DEG have the same approximate
! degree, namely, DEG = DEGREE(I)). If the list DEG is empty then
! HEAD(DEG) = 0.
! During supervariable detection HEAD(HASH) also serves as a
! pointer to a hash bucket.
! If HEAD(HASH) .GT. 0, there is a degree list of degree HASH. The
! hash bucket head pointer is LAST(HEAD(HASH)).
! If HEAD(HASH) = 0, then the degree list and hash bucket are
! both empty.
! If HEAD(HASH) .LT. 0, then the degree list is empty, and
! -HEAD(HASH) is the head of the hash bucket.
! After supervariable detection is complete, all hash buckets are
! empty, and the (LAST(HEAD(HASH)) = 0) condition is restored for
! the non-empty degree lists.

! DENXT(I)  For supervariable I, DENXT(I) is
! the supervariable following I in a link list, or zero if I is
! the last in the list. Used for two kinds of lists: degree lists
! and hash buckets (a supervariable can be in only one kind of
! list at a time). For element E, DENXT(E) is the number of
! variables with dense or full rows in the element E.

! W(I) The flag array W determines the status
! of elements and variables, and the external degree of elements.
! For elements:
! if W(E) = 0, then the element E is absorbed.
! if W(E) .GE. WFLG, then W(E)-WFLG is the size of the set
! |Le \ Lme|, in terms of nonzeros (the sum of ABS(NV(I))
! for each principal variable I that is both in the
! pattern of element E and NOT in the pattern of the
! current pivot element, ME).
! if WFLG .GT. WE(E) .GT. 0, then E is not absorbed and has
! not yet been seen in the scan of the element lists in
! the computation of |Le\Lme| in loop 150 below.
! ***SD: change comment to remove reference to label***
! For variables:
! during supervariable detection, if W(J) .NE. WFLG then J is
! not in the pattern of variable I.
! The W array is initialized by setting W(I) = 1 for all I, and by
! setting WFLG = 2. It is reinitialized if WFLG becomes too large
! (to ensure that WFLG+N does not cause integer overflow).


! Local variables:
! ---------------

! DEG:        the degree of a variable or element
! DEGME:      size (no. of variables), |Lme|, of the current element,
! ME (= DEGREE(ME))
! DEXT:       external degree, |Le \ Lme|, of some element E
! DMAX:       largest |Le| seen so far
! E:          an element
! ELENME:     the length, ELEN(ME), of element list of pivotal var.
! ELN:        the length, ELEN(...), of an element list
! EMP1:       stores number of rows found to be empty
! HASH:       the computed value of the hash function
! HMOD:       the hash function is computed modulo HMOD = MAX(1,N-1)
! I:          a supervariable
! IDUMMY:     loop counter
! ILAST:      the entry in a link list preceding I
! INEXT:      the entry in a link list following I
! IOVFLO:     local copy of ICNTL(5)
! J:          a supervariable
! JDUMMY:     loop counter
! JLAST:      the entry in a link list preceding J
! JNEXT:      the entry in a link list, or path, following J
! K:          the pivot order of an element or variable
! KNT1:       loop counter used during element construction
! KNT2:       loop counter used during element construction
! KNT3:       loop counter used during element construction
! LASTD:      index of the last row in the list of dense rows
! LENJ:       LEN(J)
! LN:         length of a supervariable list
! MAXMEM:     amount of memory needed for no compressions
! ME:         current supervariable being eliminated, and the
! current element created by eliminating that
! supervariable
! MEM:        memory in use assuming no compressions have occurred
! MINDEG:     current approximate minimum degree
! MU:         current mean of external degrees during dense detection
! NBD:        total number of dense rows selected
! NCMPA:      counter for the number of times IW was compressed
! NEL:        number of pivots selected so far
! NEWMEM:     amount of new memory needed for current pivot element
! NLEFT:      N-NEL, the number of nonpivotal rows/columns remaining
! NRLADU:     counter for the forecast number of reals in matrix factor
! NVI:        the number of variables in a supervariable I (= NV(I))
! NVJ:        the number of variables in a supervariable J (= NV(J))
! NVPIV:      number of pivots in current element
! P:          pointer into lots of things
! P1:         pe (i) for some variable i (start of element list)
! P2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
! P3:         index of first supervariable in clean list
! PJ:         pointer into an element or variable
! PDST:       destination pointer, for compression
! PEND:       end of memory to compress
! PME:        pointer into the current element (PME1...PME2)
! PME1:       the current element, ME, is stored in IW(PME1...PME2)
! PME2:       the end of the current element
! PN:         pointer into a "clean" variable, also used to compress
! PSRC:       source pointer, for compression
! SDEN:       used to remember whether dense rows occur
! SLENME:     number of variables in variable list of pivotal variable
! THRESH:     local copy of ICNTL(4)
! THRESM :    local integer holding the threshold used to detect quasi
! dense rows. When quasi dense rows are reintegrated in the
! graph to be processed then THRESM is modified.
! WE:         W(E)
! WFLG:       used for flagging the W array.  See description of W.
! WNVI:       WFLG-NV(I)
! X:          either a supervariable or an element

! OPS:        counter for forecast number of flops
! RELDEN :    holds average density to set THRESM automatically

! IDENSE is true if supervariable I is dense

! -------------------------------------------------------------------
! FUNCTIONS CALLED:
! -------------------------------------------------------------------

! ====================================================================
! INITIALIZATIONS
! ====================================================================

! ..
! .. Local Arrays ..
        INTEGER, ALLOCATABLE, DIMENSION (:) :: nv, last, degree, head, denxt, &
          w

! ..
! .. Local Scalars ..
        REAL (myreal_mc68) mu, relden
        INTEGER deg, degme, dext, dmax, e, elenme, eln, emp1, hash, hmod, i, &
          idummy, ilast, inext, iovflo, j, jdummy, jlast, jnext, k, knt1, &
          knt2, knt3, lastd, lenj, ln, maxmem, me, mem, mindeg, nbd, ncmpa, &
          ndme, nel, newmem, nleft, nvi, nvj, nvpiv, p, p1, p2, p3, pdst, &
          pend, pj, pme, pme1, pme2, pn, psrc, sden, slenme, thresh, thresm, &
          we, wflg, wnvi, x, l
        LOGICAL idense
! ..
! .. Intrinsic Functions ..
        INTRINSIC abs, int, log, max, min, mod, real, sqrt
! ..

        dmax = 0
        hmod = max(1,n-1)
        iovflo = icntl(5)
        lastd = 0
        mem = pfree - 1
        maxmem = mem
        mindeg = 1
        nbd = 0
        ncmpa = 0
        nel = 0
        thresh = icntl(4)
        wflg = 2
        emp1 = 0
        sden = 0

        ALLOCATE (nv(n),last(n),degree(n),head(n),denxt(n),w(n), &
          STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_alloc
          RETURN
        END IF

        IF (thresh>0) THEN
          thresm = int(sqrt(real(n)))

! ---------------------------------------------------------------------
! initialize head(n) and next(n)
! ---------------------------------------------------------------------
          relden = 0
          head(1:n) = 0
          denxt(1:n) = 0

! ---------------------------------------------------------------------
! create the degree hash buckets and linked lists
! for the dense nodes
! ---------------------------------------------------------------------
          relden = 0
          degree(1:n) = len(1:n)
          DO i = 1, n
            deg = len(i)
            IF (deg>0) THEN
              sden = 1
              relden = relden + deg
! insert node in degree list
              denxt(i) = head(deg)
              head(deg) = i
            ELSE
              emp1 = emp1 + 1
            END IF
          END DO

          IF (n==emp1) THEN
            mu = 0
          ELSE
            mu = relden/(n-emp1)
          END IF

! ---------------------------------------------------------------------

! 1) Recalculate the degree length of all nodes adjacent to
! the dense nodes in the degree list.  (Note:  Many of the
! dense nodes in the degree list will no longer be dense after
! this section.)

! 2) Constuct the ordering for the nodes not sent to AMD by
! selecting the most dense node in the degree list and
! then reduce the lengths of all adjacent nodes. Repeat this
! until no nodes are left with length higher than dense.
! The dense nodes are placed in the last(n) array.
! NOTE:  1) nodes are placed after the final value
! of lastnode in the last(n) array
! 2) the AMD routine will not effect anything after
! last node in the last(n) array.
! 3) nodes are saved in degree order and in their
! original state, i.e., no reverse mapping is
! needed on these.
! ---------------------------------------------------------------------

          IF (sden==1) THEN
            sden = 0
            k = n
10          CONTINUE

! ** get node from bucket
            me = head(k)

! ** main loop control
            IF (mu/=zero .AND. (n-1/=sden+emp1) .AND. (n-2/=sden+emp1)) THEN
              IF (me==0) THEN
                k = k - 1
                IF ((mu-((relden-2*k)/real(n-sden-1-emp1))>=(40.0_wp* &
                  log(real(n-sden-1-emp1)))/real(n-sden-1-emp1)+mu/real(n- &
                  sden-2-emp1)) .AND. k/=0) GO TO 10
              ELSE

! ** remove node from bucket
                head(k) = denxt(me)

! ** get degree of current node
                deg = degree(me)

! ** skip this node if degree was changed to less than dense
                IF (deg>=1) THEN

! ** check if degree was changed
                  IF (deg<k .AND. deg>0) THEN

! ** insert back into linked list at the lower degree
                    denxt(me) = head(deg)
                    head(deg) = me
                  ELSE

! ** update degree lengths of adjacent nodes
                    p1 = pe(me) + len(me) - 1
                    DO i = pe(me), p1
                      j = iw(i)
                      IF (degree(j)==1) THEN
                        degree(j) = 0
                        emp1 = emp1 + 1
                      ELSE IF (degree(j)<2*n+1) THEN
                        degree(j) = degree(j) - 1
                      END IF
                    END DO
                    sden = sden + 1
                    relden = relden - 2*deg
                    IF (n==sden+emp1) THEN
                      mu = 0
                    ELSE
                      mu = relden/(n-sden-emp1)
                    END IF
! Mark as dense row
                    degree(me) = deg + 2*n + 1
                    pe(me) = 0
                  END IF
                END IF
                GO TO 10
              END IF
            END IF
          END IF

          DO i = 1, n
            IF (degree(i)<2*n+1) THEN
              IF (degree(i)>0) THEN
                p1 = pe(i)
                p2 = pe(i) + len(i) - 1
                pj = p1
                DO j = p1, p2
! Search I for dense variables
                  IF (degree(iw(j))>2*n+1) THEN
                    len(i) = len(i) - 1
                  ELSE
                    iw(pj) = iw(j)
                    pj = pj + 1
                  END IF
                END DO
                degree(i) = len(i)
                IF (len(i)==0) pe(i) = 0
                len(i) = pj - p1
              ELSE
                len(i) = 0
                pe(i) = 0
              END IF
            END IF
          END DO
        END IF

        IF (thresh>0) THEN
! ----------------------------------------------------------
! initialize arrays and eliminate rows with no off-diag. nz.
! ----------------------------------------------------------
          last(1:n) = 0
          head(1:n) = 0
          denxt(1:n) = 0
          nv(1:n) = 1
          DO i = 1, n
            IF (len(i)==0 .AND. degree(i)==0) THEN
              nel = nel + 1
              elen(i) = -nel
              pe(i) = 0
              w(i) = 0
            ELSE
              w(i) = 1
              elen(i) = 0
            END IF
          END DO

          IF (n==nel) THEN
            GO TO 100
          ELSE
            thresm = n + 1
          END IF
        ELSE

          thresm = thresh
          last(1:n) = 0
          head(1:n) = 0
          nv(1:n) = 1
          degree(1:n) = len(1:n)
          DO i = 1, n
            IF (degree(i)==0) THEN
              nel = nel + 1
              elen(i) = -nel
              pe(i) = 0
              w(i) = 0
            ELSE
              w(i) = 1
              elen(i) = 0
            END IF
          END DO
        END IF

! ----------------------------------------------------------------
! initialize degree lists
! ----------------------------------------------------------------
        DO i = 1, n
          deg = degree(i)
          IF (deg>0) THEN
! ----------------------------------------------------------
! place i in the degree list corresponding to its degree
! or in the dense row list if i is dense
! ----------------------------------------------------------
! test for row density
            IF ((thresm>=0) .AND. (deg>=2*n+1)) THEN
! I is dense and will be inserted in the degree
! list of N
              deg = n
              inext = head(deg)
              IF (inext/=0) last(inext) = i
              denxt(i) = inext
              head(deg) = i
              last(i) = 0
              IF (lastd==0) lastd = i
            ELSE
! place i in the degree list corresponding to its degree
              inext = head(deg)
              IF (inext/=0) last(inext) = i
              denxt(i) = inext
              head(deg) = i
            END IF
          END IF
        END DO
        nbd = sden

! We suppress dense row selection if none of them was found in A
! in the 1st pass
        IF (nbd==0 .AND. thresh>0) thresm = -1

        DO WHILE (nel<n)

! ==================================================================
! GET PIVOT OF MINIMUM APPROXIMATE DEGREE
! ==================================================================
! -------------------------------------------------------------
! find next supervariable for elimination
! -------------------------------------------------------------
          DO deg = mindeg, n
            me = head(deg)
            IF (me>0) GO TO 20
          END DO
20        mindeg = deg

          IF (deg<n) THEN
! -------------------------------------------------------------
! remove chosen variable from linked list
! -------------------------------------------------------------
            inext = denxt(me)
            IF (inext/=0) last(inext) = 0
            head(deg) = inext
! -------------------------------------------------------------
! me represents the elimination of pivots nel+1 to nel+nv(me).
! place me itself as the first in this set.  It will be moved
! to the nel+nv(me) position when the permutation vectors are
! computed.
! -------------------------------------------------------------
            elenme = elen(me)
            elen(me) = -(nel+1)
            nvpiv = nv(me)
            nel = nel + nvpiv
            denxt(me) = 0

! ====================================================================
! CONSTRUCT NEW ELEMENT
! ====================================================================

! -------------------------------------------------------------
! At this point, me is the pivotal supervariable.  It will be
! converted into the current element.  Scan list of the
! pivotal supervariable, me, setting tree pointers and
! constructing new list of supervariables for the new element,
! me.  p is a pointer to the current position in the old list.
! -------------------------------------------------------------

! flag the variable "me" as being in the front by negating nv(me)
            nv(me) = -nvpiv
            degme = 0
            IF (elenme==0) THEN
! ----------------------------------------------------------
! There are no elements involved.
! Construct the new element in place.
! ----------------------------------------------------------
              pme1 = pe(me)
              pme2 = pme1 - 1
              DO p = pme1, pme1 + len(me) - 1
                i = iw(p)
                nvi = nv(i)
                IF (nvi>0) THEN
! ----------------------------------------------------
! i is a principal variable not yet placed in the
! generated element. Store i in new list
! ----------------------------------------------------
                  degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                  nv(i) = -nvi
                  pme2 = pme2 + 1
                  iw(pme2) = i

! ----------------------------------------------------
! remove variable i from degree list.
! ----------------------------------------------------
                  ilast = last(i)
                  inext = denxt(i)
                  IF (inext/=0) last(inext) = ilast
                  IF (ilast/=0) THEN
                    denxt(ilast) = inext
                  ELSE
! i is at the head of the degree list
                    head(degree(i)) = inext
                  END IF
                END IF
              END DO
! this element takes no new memory in iw:
              newmem = 0
            ELSE
! ----------------------------------------------------------
! construct the new element in empty space, iw (pfree ...)
! ----------------------------------------------------------
              p = pe(me)
              pme1 = pfree
              slenme = len(me) - elenme
              DO knt1 = 1, elenme
! search the elements in me.
                e = iw(p)
                p = p + 1
                pj = pe(e)
                ln = len(e)
! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
                DO knt2 = 1, ln
                  i = iw(pj)
                  pj = pj + 1
                  nvi = nv(i)
                  IF (nvi>0) THEN
! -------------------------------------------------
! compress iw, if necessary
! -------------------------------------------------
                    IF (pfree>iwlen) THEN
! prepare for compressing iw by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
! ***** SD: Seperate compression subroutine tried
! but found to be inefficient in comparison ****
                      pe(me) = p
                      len(me) = len(me) - knt1
! Check if anything left in supervariable ME
                      IF (len(me)==0) pe(me) = 0
                      pe(e) = pj
                      len(e) = ln - knt2
! Check if anything left in element E
                      IF (len(e)==0) pe(e) = 0
                      ncmpa = ncmpa + 1
! store first item in pe
! set first entry to -item
                      DO j = 1, n
                        pn = pe(j)
                        IF (pn>0) THEN
                          pe(j) = iw(pn)
                          iw(pn) = -j
                        END IF
                      END DO

! psrc/pdst point to source/destination
                      pdst = 1
                      psrc = 1
                      pend = pme1 - 1

! while loop:
                      DO idummy = 1, iwlen
                        IF (psrc>pend) THEN
                          GO TO 30
                        ELSE
! search for next negative entry
                          j = -iw(psrc)
                          psrc = psrc + 1
                          IF (j>0) THEN
                            iw(pdst) = pe(j)
                            pe(j) = pdst
                            pdst = pdst + 1
! copy from source to destination
                            lenj = len(j)
                            DO knt3 = 0, lenj - 2
                              iw(pdst+knt3) = iw(psrc+knt3)
                            END DO
                            pdst = pdst + lenj - 1
                            psrc = psrc + lenj - 1
                          END IF
                        END IF
                      END DO

! move the new partially-constructed element
30                    p1 = pdst
                      DO psrc = pme1, pfree - 1
                        iw(pdst) = iw(psrc)
                        pdst = pdst + 1
                      END DO
                      pme1 = p1
                      pfree = pdst
                      pj = pe(e)
                      p = pe(me)
                    END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                    degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                    nv(i) = -nvi
                    iw(pfree) = i
                    pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                    ilast = last(i)
                    inext = denxt(i)
                    IF (inext/=0) last(inext) = ilast
                    IF (ilast/=0) THEN
                      denxt(ilast) = inext
                    ELSE
! i is at the head of the degree list
                      head(degree(i)) = inext
                    END IF
                  END IF
                END DO

! set tree pointer and flag to indicate element e is
! absorbed into new element me (the parent of e is me)
                pe(e) = -me
                w(e) = 0
              END DO

! search the supervariables in me.
              knt1 = elenme + 1
              e = me
              pj = p
              ln = slenme

! -------------------------------------------------------
! search for different supervariables and add them to the
! new list, compressing when necessary.
! -------------------------------------------------------
              DO knt2 = 1, ln
                i = iw(pj)
                pj = pj + 1
                nvi = nv(i)
                IF (nvi>0) THEN
! -------------------------------------------------
! compress iw, if necessary
! -------------------------------------------------
                  IF (pfree>iwlen) THEN
! prepare for compressing iw by adjusting
! pointers and lengths so that the lists being
! searched in the inner and outer loops contain
! only the remaining entries.
                    pe(me) = p
                    len(me) = len(me) - knt1
! Check if anything left in supervariable ME
                    IF (len(me)==0) pe(me) = 0
                    pe(e) = pj
                    len(e) = ln - knt2
! Check if anything left in element E
                    IF (len(e)==0) pe(e) = 0
                    ncmpa = ncmpa + 1
! store first item in pe
! set first entry to -item
                    DO j = 1, n
                      pn = pe(j)
                      IF (pn>0) THEN
                        pe(j) = iw(pn)
                        iw(pn) = -j
                      END IF
                    END DO

! psrc/pdst point to source/destination
                    pdst = 1
                    psrc = 1
                    pend = pme1 - 1

! while loop:
! 122              CONTINUE
                    DO idummy = 1, iwlen
                      IF (psrc>pend) THEN
                        GO TO 40
                      ELSE
! search for next negative entry
                        j = -iw(psrc)
                        psrc = psrc + 1
                        IF (j>0) THEN
                          iw(pdst) = pe(j)
                          pe(j) = pdst
                          pdst = pdst + 1
! copy from source to destination
                          lenj = len(j)
                          DO knt3 = 0, lenj - 2
                            iw(pdst+knt3) = iw(psrc+knt3)
                          END DO
                          pdst = pdst + lenj - 1
                          psrc = psrc + lenj - 1
                        END IF
                      END IF
                    END DO

! move the new partially-constructed element
40                  p1 = pdst
                    DO psrc = pme1, pfree - 1
                      iw(pdst) = iw(psrc)
                      pdst = pdst + 1
                    END DO
                    pme1 = p1
                    pfree = pdst
                    pj = pe(e)
                    p = pe(me)
                  END IF

! -------------------------------------------------
! i is a principal variable not yet placed in Lme
! store i in new list
! -------------------------------------------------
                  degme = degme + nvi
! flag i as being in Lme by negating nv (i)
                  nv(i) = -nvi
                  iw(pfree) = i
                  pfree = pfree + 1

! -------------------------------------------------
! remove variable i from degree link list
! -------------------------------------------------
                  ilast = last(i)
                  inext = denxt(i)
                  IF (inext/=0) last(inext) = ilast
                  IF (ilast/=0) THEN
                    denxt(ilast) = inext
                  ELSE
! i is at the head of the degree list
                    head(degree(i)) = inext
                  END IF
                END IF
              END DO

              pme2 = pfree - 1
! this element takes newmem new memory in iw (possibly zero)
              newmem = pfree - pme1
              mem = mem + newmem
              maxmem = max(maxmem,mem)
            END IF

! -------------------------------------------------------------
! me has now been converted into an element in iw (pme1..pme2)
! -------------------------------------------------------------
! degme holds the external degree of new element
            degree(me) = degme
            pe(me) = pme1
            len(me) = pme2 - pme1 + 1

! -------------------------------------------------------------
! make sure that wflg is not too large.  With the current
! value of wflg, wflg+n must not cause integer overflow
! -------------------------------------------------------------
            IF (wflg>iovflo-n) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF

! ====================================================================
! COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
! where G' is the subgraph of G containing just the sparse rows)
! ====================================================================
! -------------------------------------------------------------
! Scan 1:  compute the external degrees of elements touched
! with respect to the current element.  That is:
! (w (e) - wflg) = |Le \ Lme|
! for each element e involving a supervariable in Lme.
! The notation Le refers to the pattern (list of
! supervariables) of a previous element e, where e is not yet
! absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
! The notation Lme refers to the pattern of the current element
! (stored in iw (pme1..pme2)).
! -------------------------------------------------------------
            DO pme = pme1, pme2
              i = iw(pme)
              eln = elen(i)
              IF (eln>0) THEN
! note that nv (i) has been negated to denote i in Lme:
                nvi = -nv(i)
                wnvi = wflg - nvi
                DO p = pe(i), pe(i) + eln - 1
                  e = iw(p)
                  we = w(e)
                  IF (we>=wflg) THEN
! unabsorbed element e has been seen in this loop
                    we = we - nvi
                  ELSE IF (we/=0) THEN
! e is an unabsorbed element - this is
! the first we have seen e in all of Scan 1
                    we = degree(e) + wnvi
                  END IF
                  w(e) = we
                END DO
              END IF
            END DO

! ====================================================================
! DEGREE UPDATE AND ELEMENT ABSORPTION
! ====================================================================

! -------------------------------------------------------------
! Scan 2:  for each sparse i in Lme, sum up the external degrees
! of each Le for the elements e appearing within i, plus the
! supervariables in i.  Place i in hash list.
! -------------------------------------------------------------

            DO pme = pme1, pme2
              i = iw(pme)
! remove absorbed elements from the list for i
              p1 = pe(i)
              p2 = p1 + elen(i) - 1
              pn = p1
              hash = 0
              deg = 0

! -------------------------------------------------------
! scan the element list associated with supervariable i
! -------------------------------------------------------
              DO p = p1, p2
                e = iw(p)
! dext = | Le | - | (Le \cap Lme)\D | - DENXT(e)
                dext = w(e) - wflg
                IF (dext>0) THEN
                  deg = deg + dext
                  iw(pn) = e
                  pn = pn + 1
                  hash = hash + e
                ELSE IF (dext==0) THEN
! aggressive absorption: e is not adjacent to me, but
! |Le(G') \ Lme(G')| is 0, so absorb it into me
                  pe(e) = -me
                  w(e) = 0
                END IF
              END DO

! count the number of elements in i (including me):
              elen(i) = pn - p1 + 1

! ----------------------------------------------------------
! scan the supervariables in the list associated with i
! ----------------------------------------------------------
              p3 = pn
              DO p = p2 + 1, p1 + len(i) - 1
                j = iw(p)
                nvj = nv(j)
                IF (nvj>0) THEN
! j is unabsorbed, and not in Lme.
! add to degree and add to new list
                  deg = deg + nvj
                  iw(pn) = j
                  pn = pn + 1
                  hash = hash + j
                END IF
              END DO

! ----------------------------------------------------------
! update the degree and check for mass elimination
! ----------------------------------------------------------
              IF (deg==0) THEN
! -------------------------------------------------------
! mass elimination - supervariable i can be eliminated
! -------------------------------------------------------
                pe(i) = -me
                nvi = -nv(i)
                degme = degme - nvi
                nvpiv = nvpiv + nvi
                nel = nel + nvi
                nv(i) = 0
                elen(i) = 0
              ELSE
! -------------------------------------------------------
! update the upper-bound degree of i
! A bound for the new external degree is the old bound plus
! the size of the generated element
! -------------------------------------------------------

! the following degree does not yet include the size
! of the current element, which is added later:
                degree(i) = min(deg,degree(i))

! -------------------------------------------------------
! add me to the list for i
! -------------------------------------------------------
! move first supervariable to end of list
                iw(pn) = iw(p3)
! move first element to end of element part of list
                iw(p3) = iw(p1)
! add new element to front of list.
                iw(p1) = me
! store the new length of the list in len (i)
                len(i) = pn - p1 + 1

! -------------------------------------------------------
! place in hash bucket.  Save hash key of i in last (i).
! -------------------------------------------------------
                hash = abs(mod(hash,hmod)) + 1
                j = head(hash)
                IF (j<=0) THEN
! the degree list is empty, hash head is -j
                  denxt(i) = -j
                  head(hash) = -i
                ELSE
! degree list is not empty - has j as its head
! last is hash head
                  denxt(i) = last(j)
                  last(j) = i
                END IF
                last(i) = hash
              END IF
            END DO
            degree(me) = degme

! -------------------------------------------------------------
! Clear the counter array, w (...), by incrementing wflg.
! -------------------------------------------------------------
            dmax = max(dmax,degme)
            wflg = wflg + dmax

! make sure that wflg+n does not cause integer overflow
            IF (wflg>=iovflo-n) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF
! at this point, w (1..n) .lt. wflg holds

! ====================================================================
! SUPERVARIABLE DETECTION
! ====================================================================
            DO pme = pme1, pme2
              i = iw(pme)
              IF ((nv(i)<0) .AND. (degree(i)<=n)) THEN
! only done for sparse rows
! replace i by head of its hash bucket, and set the hash
! bucket header to zero

! -------------------------------------------------------
! examine all hash buckets with 2 or more variables.  We
! do this by examing all unique hash keys for super-
! variables in the pattern Lme of the current element, me
! -------------------------------------------------------
                hash = last(i)
! let i = head of hash bucket, and empty the hash bucket
                j = head(hash)
                IF (j/=0) THEN
                  IF (j<0) THEN
! degree list is empty
                    i = -j
                    head(hash) = 0
                  ELSE
! degree list is not empty, restore last () of head
                    i = last(j)
                    last(j) = 0
                  END IF
                  IF (i/=0) THEN

! while loop:
                    DO jdummy = 1, n
                      IF (denxt(i)==0) THEN
                        GO TO 80
                      ELSE
! ----------------------------------------------------
! this bucket has one or more variables following i.
! scan all of them to see if i can absorb any entries
! that follow i in hash bucket.  Scatter i into w.
! ----------------------------------------------------
                        ln = len(i)
                        eln = elen(i)
! do not flag the first element in the list (me)
                        DO p = pe(i) + 1, pe(i) + ln - 1
                          w(iw(p)) = wflg
                        END DO

! ----------------------------------------------------
! scan every other entry j following i in bucket
! ----------------------------------------------------
                        jlast = i
                        j = denxt(i)

! while loop:
                        DO idummy = 1, n
                          IF (j==0) THEN
                            GO TO 70
                          ELSE

! -------------------------------------------------
! check if j and i have identical nonzero pattern
! -------------------------------------------------
! jump if i and j do not have same size data structure
                            IF (len(j)==ln) THEN
! jump if i and j do not have same number adj elts
                              IF (elen(j)==eln) THEN
! do not flag the first element in the list (me)

                                DO p = pe(j) + 1, pe(j) + ln - 1
! jump if an entry (iw(p)) is in j but not in i
                                  IF (w(iw(p))/=wflg) GO TO 50
                                END DO

! -------------------------------------------------
! found it!  j can be absorbed into i
! -------------------------------------------------
                                pe(j) = -i
! both nv (i) and nv (j) are negated since they
! are in Lme, and the absolute values of each
! are the number of variables in i and j:
                                nv(i) = nv(i) + nv(j)
                                nv(j) = 0
                                elen(j) = 0
! delete j from hash bucket
                                j = denxt(j)
                                denxt(jlast) = j
                                GO TO 60
                              END IF
                            END IF

! -------------------------------------------------
50                          CONTINUE
! j cannot be absorbed into i
! -------------------------------------------------
                            jlast = j
                            j = denxt(j)
                          END IF
60                        CONTINUE
                        END DO

! ----------------------------------------------------
! no more variables can be absorbed into i
! go to next i in bucket and clear flag array
! ----------------------------------------------------
70                      wflg = wflg + 1
                        i = denxt(i)
                        IF (i==0) GO TO 80
                      END IF
                    END DO
                  END IF
                END IF
              END IF
80            CONTINUE
            END DO

! ====================================================================
! RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
! Squeeze out absorbed variables
! ====================================================================
            p = pme1
            nleft = n - nel
            DO pme = pme1, pme2
              i = iw(pme)
              nvi = -nv(i)
              IF (nvi>0) THEN
! i is a principal variable in Lme
! restore nv (i) to signify that i is principal
                nv(i) = nvi
                IF (degree(i)<=n) THEN
! -------------------------------------------------------
! compute the external degree (add size of current elem)
! -------------------------------------------------------
                  deg = min(degree(i)+degme-nvi,nleft-nvi)
                  degree(i) = deg
                  idense = .FALSE.

! -------------------------------------------------------
! place the supervariable at the head of the degree list
! -------------------------------------------------------
                  inext = head(deg)
                  IF (inext/=0) last(inext) = i
                  denxt(i) = inext
                  last(i) = 0
                  head(deg) = i
! -------------------------------------------------------
! save the new degree, and find the minimum degree
! -------------------------------------------------------
                  mindeg = min(mindeg,deg)
                END IF
! -------------------------------------------------------
! place the supervariable in the element pattern
! -------------------------------------------------------
                iw(p) = i
                p = p + 1
              END IF
            END DO

! =====================================================================
! FINALIZE THE NEW ELEMENT
! =====================================================================
            nv(me) = nvpiv + degme
! nv (me) is now the degree of pivot (including diagonal part)
! save the length of the list for the new element me
            len(me) = p - pme1
            IF (len(me)==0) THEN
! there is nothing left of the current pivot element
              pe(me) = 0
              w(me) = 0
            END IF
            IF (newmem/=0) THEN
! element was not constructed in place: deallocate part
! of it (final size is less than or equal to newmem,
! since newly nonprincipal variables have been removed).
              pfree = p
              mem = mem - newmem + len(me)
            END IF

! =====================================================================
! END WHILE (selecting pivots)
          ELSE
! DEGREE(ME).GT.N+1 so ME is dense
! RESTARTING STRATEGY
! FOR EACH  dense row d
! 1/ insert d in the degree list according to the
! value degree(d)-(N+1) (updating MINDEG)
! 2/ ME is assumed to have no adjacent variables because just
! sorting according to order removed from matrix at initialisation
! 3/ get back to min degree process

! While loop: ME is the current dense row
! make sure that WFLG is not too large
            IF (wflg>iovflo-nbd-1) THEN
              DO x = 1, n
                IF (w(x)/=0) w(x) = 1
              END DO
              wflg = 2
            END IF
            wflg = wflg + 1
            DO idummy = 1, n
              pe(me) = 0
              len(me) = 0

! ---------------------------------------------------------
! remove chosen variable from link list
! ---------------------------------------------------------
              inext = denxt(me)
              IF (inext/=0) THEN
                last(inext) = 0
              ELSE
                lastd = 0
              END IF
! ----------------------------------------------------------
! build adjacency list of ME in quotient graph
! and calculate its external degree in ndense(me)
! ----------------------------------------------------------
              denxt(me) = 0
! Flag ME as having been considered in this calculation
              w(me) = wflg
              p1 = pe(me)
              p2 = p1 + len(me) - 1
! LN-1 holds the pointer in IW to last elt/var in adj list
! of ME.  LEN(ME) will then be set to LN-P1
! ELN-1 hold the pointer in IW to  last elt in in adj list
! of ME.  ELEN(ME) will then be set to ELN-P1
! element adjacent to ME
              ln = p1
              eln = p1

! ----------------------------------------------
! DEGREE(ME)-(2*N+1) holds last external degree computed
! when ME was detected as dense
! DENXT(ME) is the exact external degree of ME
! ----------------------------------------------
              wflg = wflg + 1
              len(me) = ln - p1
              elen(me) = eln - p1
              ndme = denxt(me) + nv(me)
              IF (denxt(me)==0) denxt(me) = 1
! ---------------------------------------------------------
! place ME in the degree list of DENXT(ME), update DEGREE
! ---------------------------------------------------------
              deg = degree(me) - 2*n - 1
              degree(me) = denxt(me)
              mindeg = min(deg,mindeg)
              jnext = head(deg)
              IF (jnext/=0) last(jnext) = me
              denxt(me) = jnext
              head(deg) = me
!  nel = nel+1

! ------------------------------
! process dense row
! ------------------------------
              me = inext
              IF (me==0) THEN
                GO TO 90
              ELSE IF (degree(me)<=(n+1)) THEN
                GO TO 90
              END IF
            END DO
90          head(n) = me
! get back to min degree elimination loop
          END IF
        END DO
! =====================================================================

100     CONTINUE
! ===================================================================
! COMPUTE THE PERMUTATION VECTORS
! ===================================================================

! ----------------------------------------------------------------
! The time taken by the following code is O(n).  At this
! point, elen (e) = -k has been done for all elements e,
! and elen (i) = 0 has been done for all nonprincipal
! variables i.  At this point, there are no principal
! supervariables left, and all elements are absorbed.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! compute the ordering of unordered nonprincipal variables
! ----------------------------------------------------------------

        l = n
        DO i = 1, n
          IF (elen(i)==0) THEN
! ----------------------------------------------------------
! i is an un-ordered row.  Traverse the tree from i until
! reaching an element, e.  The element, e, was the
! principal supervariable of i and all nodes in the path
! from i to when e was selected as pivot.
! ----------------------------------------------------------
            j = -pe(i)
! while (j is a variable) do:
            DO jdummy = 1, n
              IF (elen(j)<0) THEN
                GO TO 110
              ELSE
                j = -pe(j)
              END IF
            END DO
110         e = j
! ----------------------------------------------------------
! get the current pivot ordering of e
! ----------------------------------------------------------
            k = -elen(e)

! ----------------------------------------------------------
! traverse the path again from i to e, and compress the
! path (all nodes point to e).  Path compression allows
! this code to compute in O(n) time.  Order the unordered
! nodes in the path, and place the element e at the end.
! ----------------------------------------------------------
            j = i
! while (j is a variable) do:
            DO idummy = 1, n
              IF (elen(j)<0) THEN
                GO TO 120
              ELSE
                jnext = -pe(j)
                pe(j) = -e
                IF (elen(j)==0) THEN
! j is an unordered row
                  elen(j) = k
                  k = k + 1
                END IF
                j = jnext
              END IF
            END DO
! leave elen (e) negative, so we know it is an element
120         elen(e) = -k
          END IF
        END DO

! ----------------------------------------------------------------
! reset the inverse permutation (elen (1..n)) to be positive,
! and compute the permutation (last (1..n)).
! ----------------------------------------------------------------
        DO i = 1, n
          k = abs(elen(i))
          last(k) = i
          elen(i) = k
        END DO

! ====================================================================
! RETURN THE MEMORY USAGE IN IW AND SET INFORMATION ARRAYS
! ====================================================================
! If maxmem is less than or equal to iwlen, then no compressions
! occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
! compressions did occur, and iwlen would have had to have been
! greater than or equal to maxmem for no compressions to occur.
! Return the value of maxmem in the pfree argument.

        DEALLOCATE (nv,last,degree,head,denxt,w,STAT=info%stat)
        IF (info%stat>0) THEN
          info%flag = mc68_err_memory_dealloc
          RETURN
        END IF

        info%n_compressions = ncmpa
        info%n_zero_eigs = -1
        info%n_dense_rows = nbd
        pfree = maxmem

      END SUBROUTINE amdd

! -------------------------------------------------------------------


    END MODULE hsl_mc68_integer

! Additional modules to provide backwards compatibility.
! These modules are deprecated and may be removed at a later date.
    module hsl_mc68_double
      use hsl_mc68_integer
    end module hsl_mc68_double
    module hsl_mc68_single
      use hsl_mc68_integer
    end module hsl_mc68_single
! COPYRIGHT (c) 2012 Science and Technology Facilities Council
! Original date 6 June 2012, Version 1.0.0
!
! Written by: Jonathan Hogg and Jennifer Scott

! Given a sparse symmetric  matrix A, this package 
! uses a matching algorithm to compute an elimination
! order that is suitable for use with a sparse direct solver. 
! It optionally computes scaling factors.

! Note: this version uses mc64 to compute the matching
!       and follows approach of Duff and Pralet (2005).
!       It does not compute the optimal matching in the
!       structurally singular case as we found the extra work
!       required to do this did not result in a better
!       pivot order.

! To convert from double to single:
! * Change wp
! * Replace _double by _single
! * Replace call to mc64wd to mc64w

module hsl_mc80_double
   use hsl_mc34_double
   use hsl_mc68_integer
   ! Also calls mc64
   implicit none

   private
   public :: mc80_order, mc80_order_full, mc80_control, mc80_info

   integer, parameter :: wp = kind(0d0)
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: rinf = huge(0.0_wp)

   ! Error flags
   integer, parameter :: MC80_SUCCESS               = 0
   integer, parameter :: MC80_ERROR_ALLOCATION      = -1
   integer, parameter :: MC80_ERROR_A_N_OOR         = -2
   integer, parameter :: MC80_ERROR_SINGULAR        = -3
   integer, parameter :: MC80_ERROR_MC68            = -4
   integer, parameter :: MC80_ERROR_ORD_OOR         = -5
   integer, parameter :: MC80_ERROR_NO_METIS        = -6

   ! warning flags
   integer, parameter :: MC80_WARNING_SINGULAR      = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Data type for information returned by code
   !
   type mc80_info
      integer :: compress_rank = 0 ! order of compressed matrix (ordering
         ! applied to this matrix)
      integer :: flag = 0 ! Takes one of the enumerated flag values:
         ! Possible  following values:
         !    0 : successful entry (for structurally nonsingular matrix).
         !   -1 : allocation error
         !   -2 : Metis not available when required
         !   -3 : either n or ord is out of range (immediate return)
         !   -4 : unexpected error returned by hsl_mc68
         !   +1 : successful entry (for structurally singular matrix).
      integer :: flag68 = 0 ! error flag from hsl_mc68
      integer :: max_cycle = 0 ! maximum cycle length
      integer :: struct_rank = 0 ! structural rank 
      integer :: stat = 0 ! stat parameter

   end type mc80_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type mc80_control
      logical :: action = .true. ! if set to .true. and the matrix
         ! is found to be structurally singular, the code exits immediately
         ! with an error flag set.
      logical :: unmatched_scale_zero = .false. ! If true, then in the singular
         ! case, the scaling factors associated with the unmatched part of the
         ! matrix are set to zero. This will breakdown if the matched submatrix
         ! is numerically singular, but if no any factorization will likely have
         ! far fewer delayed pivots. If false, then Duff and Pralet are followed
         ! and the value is set such that the corresponding scaled entries are
         ! <= 1 in absolute value.
      logical :: unmatched_last = .false. ! If true, then in singular case, rows
         ! and columns associated with the unmatched part are ordered last in
         ! the elimination order. If false, they are placed in a fill-minimising
         ! position.
   end type mc80_control

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   interface mc80_order
      ! to be used if lower triangle of A available
      module procedure mc80_order_double
   end interface

   interface mc80_order_full
      ! to be used if lower and upper triangles of A available
      module procedure mc80_order_full_double
   end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine computes the matching, splits long cycles,
! compresses the matrix, applies AMD, Min Deg or Metis to the compressed
! matrix and then returns an ordering for the original
! matrix. This ordering flags 2x2 pivots using negative signs
! (as used on input to hsl_ma77).
! If scale is present, the mc64 scaling factors are returned.
! The matrix may be singular.
!
! Input (ptr, row , val) is the ** lower triangular part ** of the matrix.
! Diagonal entries need not be present.
! There is no matrix data checking.
!
subroutine mc80_order_double(ord, n, ptr, row, val, order, &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix.
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
      info%flag = MC80_ERROR_A_N_OOR
      return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
      info%flag = MC80_ERROR_ORD_OOR
      return
   end if

   ! just return with no action if n = 0
   if (n.eq. 0) return

   !
   ! Take absolute values, expand out, removing any explicit zeroes
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(2*ne), val2(2*ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc34_expand(n, row2, ptr2, cperm, a=val2)
   
   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
           perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
           perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! This routine is the same as mc80_order EXCEPT on
! input ptr, row , val hold the ** lower AND upper ** triangular 
! parts of the matrix.
! this reduces amount of copies of matrix required (so slightly
! more efficient on memory and does not need to expand supplied matrix)
!  
subroutine mc80_order_full_double(ord, n, ptr, row, val, order,  &
      control, info, scale)
   integer, intent(in) :: ord ! controls ordering on compressed matrix
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   integer, dimension(:), intent(out) :: order ! |order(i)|  holds the position
      ! of variable i in the elimination order (pivot sequence). If a
      ! 1x1 pivot i is obtained,  order(i)>0. If a
      ! 2x2 pivot involving  i and j is obtained, 
      ! order(i)<0, order(j)<0 and |order(j)|=|order(i)|+1. 
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(out) :: info ! used to hold information
   real(wp), dimension(n), intent(out), optional :: scale ! if present,
      ! returns the mc64 symmetric scaling

   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: ptr2 ! column pointers for expanded 
      ! matrix. 
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix 
   real(wp), dimension(:), allocatable :: val2 ! entries of expanded matrix.
   real(wp), dimension(:), allocatable :: scale2 ! holds scaling factors
      ! if scale is not present.

   integer :: i, j, k, ne

   info%compress_rank = 0
   info%flag = 0
   info%flag68 = 0
   info%stat = 0
   info%max_cycle = 0 
   info%struct_rank = n

   ! check n has valid value
   if (n < 0) then
     info%flag = MC80_ERROR_A_N_OOR
     return
   end if

   ! check ord has valid value
   if (ord < 1 .or. ord > 3) then
     info%flag = MC80_ERROR_ORD_OOR
     return
   end if

   ! just return with no action if n = 0
   if (n.eq.0) return

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne), cperm(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = abs(val(j))
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   ! Compute matching and scaling

   if (present(scale)) then 
      call mc80_scale(n,ptr2,row2,val2,scale,control,info, &
         perm=cperm)
   else
      allocate(scale2(n), stat=info%stat)
      if (info%stat.ne.0) then
         info%flag = MC80_ERROR_ALLOCATION
         return
      end if
      call mc80_scale(n,ptr2,row2,val2,scale2,control,info, &
         perm=cperm)
      deallocate(scale2, stat=info%stat)
   end if
   deallocate(val2, stat=info%stat)

   if (info%flag.lt.0) return

   ! Note: row j is matched with column cperm(j)
   ! write (*,'(a,15i4)') 'cperm',cperm(1:min(15,n))
   !
   ! Split matching into 1- and 2-cycles only and then
   ! compress matrix and order.

   call mc80_split(ord,n,row2,ptr2,order,cperm,control,info)

   if(present(scale)) scale(1:n) = exp( scale(1:n) )

end subroutine mc80_order_full_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Split matching into 1- and 2-cycles only and then
! compress matrix and order using mc68.
!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! Overwritten in the singular case
!
subroutine mc80_split(ord,n,row2,ptr2,order,cperm,control,info)
   integer, intent(in) :: ord ! controls choice of ordering algorithm 
      ! on compressed matrix
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr2 
   integer, dimension(:), intent(in) :: row2 
   integer, dimension(n), intent(out) :: order ! used to hold ordering
   integer, dimension(n), intent(inout) :: cperm ! used to hold matching 
   type (mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information

   type (mc68_control) :: control68
   type (mc68_info) :: info68

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in condensed
      ! matrix.
   integer, dimension(:), allocatable :: ptr3 ! column pointers for condensed 
      ! matrix.
   integer, dimension(:), allocatable :: row3 ! row indices for condensed 
      ! matrix.


   integer :: csz ! current cycle length
   integer :: i, j, j1, j2, jj, k, krow
   integer :: max_csz ! maximum cycle length
   integer :: ncomp ! order of compressed matrix
   integer :: ncomp_matched ! order of compressed matrix (matched entries only)
   integer :: ne ! number of non zeros

   ! Use iwork to track what has been matched:
   ! -2 unmatched
   ! -1 matched as singleton
   !  0 not yet seen
   ! >0 matched with specified node

   ne = ptr2(n+1) - 1
   allocate(ptr3(n+1), row3(ne), old_to_new(n), new_to_old(n), iwork(n), &
      stat=info%stat)
   if (info%stat.ne.0) return

   iwork(1:n) = 0
   max_csz = 0
   do i = 1, n
      if (iwork(i).ne.0) cycle
      j = i
      csz = 0
      do
         if (cperm(j).eq.-1) then
            ! unmatched by MC64
            iwork(j) = -2
            csz = csz + 1
            exit
         else if (cperm(j).eq.i) then
            ! match as singleton, unmatched or finished
            iwork(j) = -1
            csz = csz + 1
            exit
         end if
         ! match j and cperm(j)
         jj = cperm(j)
         iwork(j) = jj
         iwork(jj) = j
         csz = csz + 2
         ! move onto next start of pair
         j = cperm(jj)
         if (j.eq.i) exit
      end do
      max_csz = max(max_csz, csz)
   end do
   info%max_cycle = max_csz

   ! Overwrite cperm with new matching
   cperm(1:n) = iwork(1:n)

   !
   ! Build maps for new numbering schemes
   !
   k = 1
   do i = 1,n
      j = cperm(i)
      if (control%unmatched_last .and. j.eq.-2) cycle
      if (j<i .and. j.gt.0) cycle
      old_to_new(i) = k
      new_to_old(k) = i ! note: new_to_old only maps to first of a pair
      if (j.gt.0) old_to_new(j) = k   
      k = k + 1
   end do
   ncomp_matched = k-1
   
   ! Place unmatched columns to rear
   if(control%unmatched_last) then
      do i = 1, n
         if(cperm(i).eq.-2) then
            old_to_new(i) = k
            new_to_old(k) = i
            k = k + 1
         endif
      end do
   end if

   !
   ! Produce a condensed version of the matrix for ordering.
   ! Hold pattern using ptr3 and row3.
   !
   ptr3(1) = 1
   iwork(:) = 0 ! Use to indicate if entry is in a paired column
   ncomp = 1
   jj = 1
   do i = 1, n
      j = cperm(i)
      if (j<i .and. j.gt.0) cycle ! already seen
      if (control%unmatched_last .and. j.eq.-2) cycle ! column not participating
      do k = ptr2(i), ptr2(i+1)-1
         krow = old_to_new(row2(k))
         if (iwork(krow).eq.i) cycle ! already added to column
         if (krow>ncomp_matched) cycle ! unmatched row not participating
         row3(jj) = krow
         jj = jj + 1
         iwork(krow) = i
      end do
      if (j.gt.0) then
         ! Also check column cperm(i)
         do k = ptr2(j), ptr2(j+1)-1
            krow = old_to_new(row2(k))
            if (iwork(krow).eq.i) cycle ! already added to column
            if (krow>ncomp_matched) cycle ! unmatched row not participating
            row3(jj) = krow
            jj = jj + 1
            iwork(krow) = i
         end do
      end if
      ptr3(ncomp+1) = jj
      ncomp = ncomp + 1
   end do
   ncomp = ncomp - 1
   info%compress_rank = ncomp

   if(control%unmatched_last) ncomp = ncomp_matched

   ! store just lower triangular part for input to hsl_mc68
   ptr3(1) = 1
   jj = 1
   j1 = 1
   do i = 1, ncomp
      j2 = ptr3(i+1)
      do k = j1, j2-1
         krow = row3(k)
         if ( krow.lt.i ) cycle ! already added to column
         row3(jj) = krow
         jj = jj + 1
      end do
      ptr3(i+1) = jj
      j1 = j2
   end do

   ! reorder the compressed matrix using hsl_mc68.
   ! switch off hsl_mc68 printing
   control68%lp = -1
   control68%wp = -1
   control68%mp = -1
   control68%print_level = -1
   call mc68_order(ord,ncomp,ptr3,row3,order,control68,info68)

   if (info68%flag < 0) then
      info%flag68 = info68%flag
      select case(info68%flag)
      case(-1)
         info%flag = MC80_ERROR_ALLOCATION
         info%stat = info68%stat
      case(-5)
         info%flag = MC80_ERROR_NO_METIS
      case default
         info%flag = MC80_ERROR_MC68
      end select
      return
   end if

   do i = 1, ncomp
      j = order(i)
      iwork(j) = i
   end do

   !
   ! Translate inverse permutation in iwork back to 
   ! permutation for original variables.
   ! Set negative signs for 2x2 pivots (exploited by hsl_ma77).
   !
   k = 1
   do i = 1, ncomp
      j = new_to_old( iwork(i) )
      order(j) = k
      k = k + 1
      if (cperm(j).gt.0) then
         order(j) = -order(j)
         j = cperm(j)
         order(j) = -k
         k = k + 1
      end if
   end do

   ! Place unmatched columns last
   if(control%unmatched_last) then
      do i = ncomp+1, ncomp + (n-info%struct_rank)
         j = new_to_old( i )
         order(j) = k
         k = k + 1
      end do
   endif

end subroutine mc80_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Scale the matrix using MC64, accounting for singular matrices using the
! approach of Duff and Pralet
!
! Expects a full matrix as input
!
subroutine mc80_scale(n, ptr, row, val, scale, control, info, perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(in) :: ptr
   integer, dimension(:), intent(in) :: row
   real(wp), dimension(:), intent(in) :: val
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type(mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: ptr2 ! column pointers after 
      ! zeros removed.
   integer, dimension(:), allocatable :: row2 ! row indices after zeros
      !  removed. 
   real(wp), dimension(:), allocatable :: val2 ! matrix of absolute values
      ! (zeros removed).
   real(wp), dimension(:), allocatable :: cscale ! temporary copy of scaling
      ! factors. Only needed if A rank deficient. allocated to have size n.

   integer :: i, j, k, ne

   info%struct_rank = n

   !
   ! Remove any explicit zeroes and take absolute values
   !
   ne = ptr(n+1) - 1
   allocate(ptr2(n+1), row2(ne), val2(ne),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 1
   do i = 1, n
      ptr2(i) = k
      do j = ptr(i), ptr(i+1)-1
         if (val(j).eq.zero) cycle
         row2(k) = row(j)
         val2(k) = val(j)
         k = k + 1
      end do
   end do
   ptr2(n+1) = k

   call mc80_match(n,row2,ptr2,val2,scale,control,info,perm=perm)
   if (info%flag.lt.0) return

   if (info%struct_rank.ne.n) then
      ! structurally singular case. At this point, scaling factors
      ! for rows in corresponding to rank deficient part are set to 
      ! zero. The following is to set them according to Duff and Pralet.
      deallocate(ptr2, stat=info%stat)
      deallocate(row2, stat=info%stat)
      deallocate(val2, stat=info%stat)
      if(.not.control%unmatched_scale_zero) then
         allocate(cscale(n),stat=info%stat)
         if (info%stat.ne.0) then
            info%flag = MC80_ERROR_ALLOCATION
            return
         end if
         cscale(1:n) = scale(1:n)
         do i = 1,n
            if (cscale(i).ne.-huge(scale)) cycle
            do j = ptr(i), ptr(i+1)-1
               k = row(j)
               if (cscale(k).eq.-huge(scale)) cycle
               scale(i) = max(scale(i), val(j)+scale(k))
            end do
            if(scale(i).eq.-huge(scale)) then
               scale(i) = zero
            else
               scale(i) = -scale(i)
            endif
         end do
      end if
   end if

end subroutine mc80_scale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Input (ptr2, row2 , val2) holds the ** lower and upper triangles ** 
! of the matrix (with explicit zeros removed).
! val2 holds absolute values of matrix entries.
! Overwritten in the singular case
!
subroutine mc80_match(n,row2,ptr2,val2,scale,control,info,perm)
   integer, intent(in) :: n
   integer, dimension(:), intent(inout) :: ptr2 ! In singular case, overwritten
      ! by column pointers for non singular part of matrix.
   integer, dimension(:), intent(inout) :: row2 ! In singular case, overwritten
      ! by row indices for non singular part of matrix.
   real(wp), dimension(:), intent(inout) :: val2 ! In singular case, overwritten
      ! by entries for non singular part of matrix.
   real(wp), dimension(n), intent(out) :: scale ! returns the symmetric scaling
   type(mc80_control), intent(in) :: control
   type (mc80_info), intent(inout) :: info ! used to hold information
   integer, dimension(n), intent(out), optional :: perm ! if present, returns
      ! the matching

   integer, dimension(:), allocatable :: iwork ! work array
   integer, dimension(:), allocatable :: cperm ! used to hold matching
   integer, dimension(:), allocatable :: old_to_new, new_to_old
      ! holds mapping between original matrix indices and those in reduced
      ! non singular matrix. 
   real(wp), dimension(:), allocatable :: cmax ! (log) column maximum
   real(wp), dimension(:), allocatable :: dw ! array used by mc64

   integer :: i, j, j1, j2, jj, k
   integer :: ne ! number of non zeros
   integer :: nn ! Holds number of rows/cols in non singular part of matrix
   integer :: nne ! Only used in singular case. Holds number of non zeros
     ! in non-singular part of matrix.
   integer :: rank ! returned by mc64
   real(wp) :: colmax ! max. entry in col. of expanded matrix

   allocate(iwork(5*n), cperm(n), dw(2*n), cmax(n), stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if
    
   ! Compute column maximums    
   do i = 1,n
      colmax = max(zero,maxval(val2(ptr2(i):ptr2(i+1)-1)))
      if (colmax.ne.zero) colmax = log(colmax)
      cmax(i) = colmax
   end do

   do i = 1,n
      val2(ptr2(i):ptr2(i+1)-1) = cmax(i) - log(val2(ptr2(i):ptr2(i+1)-1))
   end do

   ne = ptr2(n+1)-1
   call mc64wd(n,ne,ptr2,row2,val2,cperm,rank, &
      iwork(1),iwork(n+1),iwork(2*n+1),iwork(3*n+1),iwork(4*n+1), &
      dw(1),dw(n+1))

   if (rank.eq.n) then
      do i = 1,n
         scale(i) = (dw(i)+dw(n+i)-cmax(i))/2
      end do
      if (present(perm)) perm(1:n) = cperm(1:n)
      return
   end if

   !!!! we have to handle the singular case. Either immediate exit
   ! or set warning, squeeze out the unmatched entries and recall mc64wd.

   info%struct_rank = rank
   if (.not.control%action) then
      info%flag = MC80_ERROR_SINGULAR
      return
   end if

   info%flag = MC80_WARNING_SINGULAR

   allocate(old_to_new(n), new_to_old(n),stat=info%stat)
   if (info%stat.ne.0) then
      info%flag = MC80_ERROR_ALLOCATION
      return
   end if

   k = 0
   do i = 1,n
      if (cperm(i) < 0) then
         ! row i and col j are not part of the matching
         old_to_new(i) = -1
      else
         k = k + 1
         ! old_to_new(i) holds the new index for variable i after
         ! removal of singular part and new_to_old(k) is the
         ! original index for k
         old_to_new(i) = k
         new_to_old(k) = i
      end if
   end do

   ! Overwrite ptr2, row2 and val2
   nne = 0
   k = 0
   ptr2(1) = 1
   j2 = 1
   do i = 1,n
      j1 = j2
      j2 = ptr2(i+1)
      ! skip over unmatched entries
      if (cperm(i) < 0) cycle
      k = k + 1
      do j = j1,j2-1
         jj = row2(j)
         if (cperm(jj) < 0) cycle
         nne = nne + 1
         row2(nne) = old_to_new(jj)
         val2(nne) = val2(j)
      end do
      ptr2(k+1) = nne + 1
    end do
    ! nn is order of non-singular part.
    nn = k
    call mc64wd(nn,nne,ptr2,row2,val2,cperm,rank, &
       iwork(1),iwork(nn+1),iwork(2*nn+1),iwork(3*nn+1),iwork(4*nn+1), &
       dw(1),dw(nn+1))

    do i = 1,n
       j = old_to_new(i)
       if (j < 0) then
          scale(i) = -huge(scale)
       else
         ! Note: we need to subtract col max using old matrix numbering
         scale(i) = (dw(j)+dw(nn+j)-cmax(i))/2
      end if
   end do

   if (present(perm)) then
      perm(1:n) = -1
      do i = 1,nn
         j = cperm(i)
         perm(new_to_old(i)) = new_to_old(j)
      end do
   end if

end subroutine mc80_match

!**********************************************************************
end module hsl_mc80_double
