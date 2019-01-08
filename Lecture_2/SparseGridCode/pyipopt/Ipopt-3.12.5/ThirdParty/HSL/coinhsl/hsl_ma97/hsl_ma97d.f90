! COPYRIGHT (c) 2011,2012 Science and Technology Facilities Council
! Original date 20 December 2011, Version 1.0.0
!
! Written by: Jonathan Hogg and Jennifer Scott
!
! Version 2.1.0
! See ChangeLog for history
!

!
! To convert from double:
! * Change wp
! * Change _double
! * Change BLAS calls: daxpy, dgemm, dgemv, dpotrf, dswap, dsyrk, dtrsm, dtrsv,
!                      dnrm2, dtrmv, dtrmm
! * Change HSL calls: mc30ad, mc77id, mc77ad
! * Change control%u default to 0.1
! * Change control%small to 1e-12
!
module hsl_ma97_double
!$ use omp_lib
   !%%use logger
   use hsl_mc34_double
   use hsl_mc64_double
   use hsl_mc68_integer
   use hsl_mc69_double
   use hsl_mc78_integer
   use hsl_mc80_double
   implicit none

   private
   public :: ma97_akeep, ma97_fkeep, ma97_control, ma97_info
   public :: ma97_analyse, ma97_analyse_coord, ma97_factor, ma97_factor_solve, &
      ma97_solve, ma97_solve_fredholm, ma97_free, ma97_finalise, &
      ma97_enquire_posdef, ma97_enquire_indef, ma97_alter, ma97_lmultiply, &
      ma97_sparse_fwd_solve
   public :: ma97_get_n__, ma97_get_nz__
   !public :: cmp_fkeep

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: zero = 0.0_wp

   integer, parameter :: nemin_default = 8 ! node amalgamation parameter

   ! Success flag
   integer, parameter :: MA97_SUCCESS               = 0

   ! Error flags
   integer, parameter :: MA97_ERROR_CALL_SEQUENCE   = -1
   integer, parameter :: MA97_ERROR_A_N_OOR         = -2
   integer, parameter :: MA97_ERROR_A_PTR           = -3
   integer, parameter :: MA97_ERROR_A_ALL_OOR       = -4
   integer, parameter :: MA97_ERROR_MATRIX_TYPE     = -5
   integer, parameter :: MA97_ERROR_IMAG_DIAGONAL   = -6
   integer, parameter :: MA97_ERROR_SINGULAR        = -7
   integer, parameter :: MA97_ERROR_NOT_POS_DEF     = -8
   integer, parameter :: MA97_ERROR_INFINITY        = -9
   integer, parameter :: MA97_ERROR_PTR_ROW         = -10
   integer, parameter :: MA97_ERROR_ORDER           = -11
   integer, parameter :: MA97_ERROR_X_SIZE          = -12
   integer, parameter :: MA97_ERROR_JOB_OOR         = -13
   integer, parameter :: MA97_ERROR_NOT_LLT         = -14
   integer, parameter :: MA97_ERROR_NOT_LDLT        = -15
   integer, parameter :: MA97_ERROR_ALLOCATION      = -16
   integer, parameter :: MA97_ERROR_NO_METIS        = -17
   integer, parameter :: MA97_ERROR_MC68            = -18
   integer, parameter :: MA97_ERROR_MC77            = -19
   integer, parameter :: MA97_ERROR_VAL             = -20
   integer, parameter :: MA97_ERROR_NO_SAVED_SCALING= -21
   integer, parameter :: MA97_ERROR_NBI_OOR         = -22
   integer, parameter :: MA97_ERROR_UNKNOWN         = -99

   ! warning flags
   integer, parameter :: MA97_WARNING_IDX_OOR          = 1
   integer, parameter :: MA97_WARNING_DUP_IDX          = 2
   integer, parameter :: MA97_WARNING_DUP_AND_OOR      = 3
   integer, parameter :: MA97_WARNING_MISSING_DIAGONAL = 4
   integer, parameter :: MA97_WARNING_MISS_DIAG_OORDUP = 5
   integer, parameter :: MA97_WARNING_ANAL_SINGULAR    = 6
   integer, parameter :: MA97_WARNING_FACT_SINGULAR    = 7
   integer, parameter :: MA97_WARNING_MATCH_ORD_NO_SCALE=8

   ! solve job values
   integer, parameter :: MA97_SOLVE_JOB_ALL        = 0 ! PLD(PL)^TX = B
   integer, parameter :: MA97_SOLVE_JOB_FWD        = 1 ! PLX = B
   integer, parameter :: MA97_SOLVE_JOB_DIAG       = 2 ! DX = B (indef only)
   integer, parameter :: MA97_SOLVE_JOB_BWD        = 3 ! (PL)^TX = B
   integer, parameter :: MA97_SOLVE_JOB_DIAG_BWD   = 4 ! D(PL)^TX = B (indef)

   ! Assorted fixed parameters (could be made part of control)
   integer(long), parameter :: BLK_SZ = 2**20 ! 8+4MB pages (real+int)

   !
   ! These parameters are used for diagnostic output for generating profiles
   ! and pictures of the parallel tree decomposition
   !
   !%%!logical, parameter :: dolog = .true.
   !%%logical, parameter :: dolog = .false.
   !%%integer, parameter :: LOG_UNIT = 15

   !%%!logical, parameter :: dograph = .true.
   !%%logical, parameter :: dograph = .false.
   !%%integer, parameter :: GRAPH_UNIT = 16

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Type for custom allocator
   ! Used to aggregate many small allocations by doing a single big allocation
   ! and chopping it up.
   ! Note: Only supports freeall operation, not individual frees.
   type smalloc_type
      real(wp), dimension(:), allocatable :: rmem ! real memory
      integer(long) :: rmem_size ! needed as size(rmem,kind=long) is f2003
      integer(long) :: rhead = 0 ! last location containing useful information
         ! in rmem
      integer, dimension(:), allocatable :: imem ! integer memory
      integer(long) :: imem_size ! needed as size(imem,kind=long) is f2003
      integer(long) :: ihead = 0 ! last location containing useful information
         ! in imem
      type(smalloc_type), pointer :: next_alloc => null()
      type(smalloc_type), pointer :: top_real => null() ! Last page where real
         ! allocation was successful
      type(smalloc_type), pointer :: top_int => null() ! Last page where integer
         ! allocation was successful
!$    integer(omp_lock_kind) :: lock
   end type smalloc_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Stack memory allocation type
   type stack_mem_type
      real(wp), dimension(:), allocatable :: mem ! real memory
      integer(long) :: mem_size ! needed as size(mem,kind=long) is f2003
      integer(long) :: head = 0 ! last location containing useful information
      type(stack_mem_type), pointer :: below => null() ! next stack frame down
   end type stack_mem_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for information generated in analyse phase
   !
   type ma97_akeep
      private

      logical :: check ! copy of check as input to analyse phase
      integer :: flag ! copy of error flag.
      integer :: maxmn ! maximum value of blkm or blkn
      integer :: n ! Dimension of matrix
      integer :: ne ! Set to number of entries input by user.
      integer(long) :: nfactor 
      integer :: nnodes = -1 ! Number of nodes in assembly tree
      integer :: num_two ! This is set to 0 as we ignore any negative signs
         ! that indicate 2x2 pivot (analyse does not exploit them)

      ! child_list(child_ptr(node):child_ptr(node+1)-1) is list of children
      ! of node. Used to ensure we always sum contributions from children
      ! in the same order
      integer, dimension(:), allocatable :: child_ptr 
      integer, dimension(:), allocatable :: child_list
      integer, dimension(:), allocatable :: invp ! inverse of pivot order that
         ! is passed to factorize phase
      integer, dimension(:), allocatable :: level ! level(i) of the assembly
         ! tree at which node i sits. root is level 1.
      integer, dimension(:,:), allocatable :: nlist ! map from A to factors
         ! For nodes i, the entries nlist(1:2, nptr(i):nptr(i+1)-1) define
         ! a relationship:
         ! nodes(node)%lcol( nlist(2,j) ) = val( nlist(1,j) )
     integer, dimension(:), allocatable :: nptr ! Entries into nlist for
         ! nodes of the assembly tree. Has length nnodes+1
      integer, dimension(:), allocatable :: rlist ! rlist(rptr(i):rptr(i+1)-1)
         ! contains the row indices for node i of the assembly tree. 
         ! At each node, the list
         ! is in elimination order. Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: rptr ! Pointers into rlist
         ! for nodes of assembly tree. Has length nnodes+1. 
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sparent ! sparent(i) is parent
         ! of node i in assembly tree. sparent(i)=nnodes+1 if i is a root.
         ! The parent is always numbered higher than each of its children.
         ! Allocated within mc78_analyse.
      integer, dimension(:), allocatable :: sptr ! (super)node pointers.
         ! Supernode i consists of sptr(i) through sptr(i+1)-1.
         ! Allocated within mc78_analyse.
      integer(long), dimension(:), allocatable :: subtree_work ! For each node,
         ! the number of flops involved in the subtree rooted at that node.

      ! Following components are for cleaned up matrix data.
      ! LOWER triangle only. We have to retain these for factorize phase
      ! as used if the user wants to do scaling.
      ! These components are NOT used if check is set to .false.
      ! on call to ma97_analyse.
      integer, allocatable :: ptr(:) ! column pointers
      integer, allocatable :: row(:)! row indices
      integer :: lmap ! used by hsl_mc69
      integer, allocatable :: map(:) ! used by hsl_mc69

      ! Following components are cached members of info
      integer :: matrix_dup
      integer :: matrix_outrange
      integer :: matrix_missing_diag
      integer :: maxdepth
      integer(long) :: num_flops ! not copied to info in factor, but used to
         ! determine if parallelism should be used
      integer :: num_sup
      integer :: ordering

      ! Scaling from matching-based ordering
      real(wp), dimension(:), allocatable :: scaling
   end type ma97_akeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Data type for storing each node of the factors
   type node_type
      integer :: nelim
      integer :: ndelay
      real(wp), dimension(:), pointer :: lcol ! values in factors
         ! (will also include unneeded data for any columns delayed from this
         ! node)
      integer, dimension(:), pointer :: perm ! permutation of columns at this
         ! node: perm(i) is column index in expected global elimination order
         ! that is actually eliminated at local elimination index i
         ! Assuming no delays or permutation this will be
         ! sptr(node):sptr(node+1)-1
      ! Following components are used to index directly into contiguous arrays
      ! lcol and perm without taking performance hit for passing pointers
      type(smalloc_type), pointer :: rsmptr, ismptr
      integer(long) :: rsmsa, ismsa
   end type node_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for data generated in factorise phase
   !
   type ma97_fkeep
      private

      integer :: flag ! copy of error flag.
      real(wp), dimension(:), allocatable :: scaling ! Stores scaling for
         ! each entry (in original matrix order)
      type(node_type), dimension(:), allocatable :: nodes ! Stores pointers
         ! to information about nodes
      type(smalloc_type), pointer :: alloc=>null() ! Linked list of memory pages
         ! pointed to by nodes variable
      logical :: pos_def ! set to true if user indicates matrix pos. definite

      ! Info components to copy on solve
      integer :: matrix_rank
      integer :: maxfront
      integer :: num_delay
      integer(long) :: num_factor
      integer(long) :: num_flops
      integer :: num_neg
      integer :: num_two
   end type ma97_fkeep

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for temporary stack data that is only needed transiently during
   ! factorise phase
   ! Each instance represents a "page" of memory
   !
   type stack_type
      real(wp), dimension(:), pointer :: val => null() ! generated element
      ! Following components allow us to pass contiguous array val without
      ! taking performance hit for passing pointers
      type(stack_mem_type), pointer :: stptr => null()
      integer(long) :: stsa
   end type stack_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for per-thread stats. This is amalgamated after end of parallel
   ! section to get info parameters of same name.
   !
   type thread_stats
      integer :: flag = MA97_SUCCESS
      integer :: st = 0
      integer :: maxfront = 0 ! Maximum front size
      integer(long) :: num_factor = 0_long ! Number of entries in factors
      integer(long) :: num_flops = 0_long ! Number of floating point operations
      integer :: num_delay = 0 ! Number of delayed variables
      integer :: num_neg = 0 ! Number of negative pivots
      integer :: num_two = 0 ! Number of 2x2 pivots
      integer :: num_zero = 0 ! Number of zero pivots
   end type thread_stats

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! This type is used to pass buf around for each thread such that it can
   ! be reallocated independantly
   !
   type real_ptr_type
      real(wp), pointer :: chkptr => null()
      real(wp), dimension(:), allocatable :: val
   end type real_ptr_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for control parameters
   !
   type ma97_control
      logical :: action = .true. ! used in indefinite case only.
         ! If true and the matrix is found to be
         ! singular, computation continues with a warning.
         ! Otherwise, terminates with error MA97_ERROR_SINGULAR.
      integer(long) :: factor_min = 20000000 ! min number of expected flops
         ! before parallel execution is used in factorization (=2e7)
      integer :: nemin = nemin_default ! Min. number of eliminations at a tree
         ! node for amalgamation not to be considered.
      real(wp) :: multiplier = 1.1 ! size to multiply expected memory size by
         ! when doing initial memory allocation to allow for delays.
      integer :: ordering = 5 ! controls choice of ordering
         ! 0 Order must be supplied by user
         ! 1 An approximate minimum degree ordering is used.
         ! 2 A minimum degree ordering is used (as in MA27).
         ! 3 METIS ordering with default settings is used.
         ! 4 MA47 ordering for indefinite matrices is used.
         ! 5 Automatic choice between AMD and METIS. Parallel  (DEFAULT)
         ! 6 Automatic choice between AMD and METIS. Serial.
         ! 7 Matching with AMD on compressed matrix.
         ! 8 Matching with METIS on compressed matrix.
      integer :: print_level = 0 ! Controls diagnostic printing.
         ! Possible values are:
         !  < 0: no printing.
         !  0: error and warning messages only.
         !  1: as 0 plus basic diagnostic printing.
         !  > 1: as 1 plus some more detailed diagnostic messages.
      integer :: scaling = 0 ! controls use of scaling. 
         !  <=0: user supplied (or no) scaling 
         !  1: mc64
         !  2: mc77 
         !  3: mc64 from mc80
         !  >=4: mc30 
      real(wp) :: small = 1d-20 ! Minimum pivot size (absolute value of a
         ! pivot must be of size at least small to be accepted).
      logical :: solve_blas3 = .false. ! Use dgemm rather than dgemv in solve
         ! with a single rhs if true
      integer(long) :: solve_min = 100000 ! Minimum value of info%num_factor to
         ! deploy parallel solve on
      logical :: solve_mf = .false. ! Do we use s/n (false) or m/f (true) solve?
      real(wp) :: u = 0.01
      integer :: unit_diagnostics = 6 ! unit number for diagnostic printing.
         ! Printing is suppressed if unit_diagnostics  <  0.
      integer :: unit_error = 6 ! unit number for error messages.
         ! Printing is suppressed if unit_error  <  0.
      integer :: unit_warning = 6 ! unit number for warning messages.
         ! Printing is suppressed if unit_warning  <  0.

      !
      ! Undocumented
      !
      integer(long) :: min_subtree_work = int(1e5) ! Minimum amount of work
         ! to aim for in each parallel task for subtree decomposition
      integer :: min_ldsrk_work = int(1e4) ! Minimum amount of work to aim
         ! for in an ldsrk subtask

      !
      ! Added after initial release
      !
      real(wp) :: consist_tol = epsilon(one) ! used on call to 
         ! ma97_solve_fredholm to
         ! determine if system is consistent (singular case only)
   end type ma97_control

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !
   ! Data type for information returned by code
   !
   type ma97_info
      integer :: flag ! Takes one of the enumerated flag values:
         ! MA97_SUCCESS
         ! MA97_ERROR_XXX
         ! MA97_WARNING_XXX
      integer :: flag68 = 0 ! error flag from hsl_mc68
      integer :: flag77 = 0 ! error flag from mc77
      integer :: matrix_dup = 0 ! Number of duplicated entries.
      integer :: matrix_rank = 0 ! Rank of matrix (anal=structral, fact=actual)
      integer :: matrix_outrange = 0 ! Number of out-of-range entries.
      integer :: matrix_missing_diag = 0 ! Number of missing diag. entries
      integer :: maxdepth ! Maximum depth of tree
      integer :: maxfront ! Maximum front size
      integer(long) :: num_factor = 0_long ! Number of entries in factors
      integer(long) :: num_flops = 0_long ! Number of floating point operations
      integer :: num_delay = 0 ! Number of delayed variables
      integer :: num_neg = 0 ! Number of negative pivots
      integer :: num_sup = 0 ! Number of supernodes
      integer :: num_two = 0 ! Number of 2x2 pivots used by factorization
      integer :: ordering = 0 ! ordering actually used
      integer :: stat = 0 ! stat parameter
   end type ma97_info

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! Make interfaces generic.
   interface ma97_analyse
      module procedure analyse_double
   end interface ma97_analyse

   interface ma97_analyse_coord
      module procedure ma97_analyse_coord_double
   end interface ma97_analyse_coord

   interface ma97_factor
      module procedure ma97_factor_double
   end interface ma97_factor

   interface ma97_factor_solve
      module procedure ma97_factor_solve_double
      module procedure ma97_factor_solve_one_double
   end interface ma97_factor_solve

   interface ma97_solve
      module procedure ma97_solve_one_double
      module procedure ma97_solve_mult_double
   end interface ma97_solve

   interface ma97_solve_fredholm
      module procedure ma97_solve_fredholm_double
   end interface ma97_solve_fredholm

   interface ma97_free
      module procedure free_akeep_double
      module procedure free_fkeep_double
   end interface ma97_free

   interface ma97_finalise
      module procedure finalise_both_double
   end interface ma97_finalise

   interface ma97_enquire_posdef
      module procedure ma97_enquire_posdef_double
   end interface ma97_enquire_posdef

   interface ma97_enquire_indef
      module procedure ma97_enquire_indef_double
   end interface ma97_enquire_indef

   interface ma97_alter
      module procedure ma97_alter_double
   end interface ma97_alter

   interface ma97_lmultiply
      module procedure ma97_lmultiply_one_double
      module procedure ma97_lmultiply_mult_double
   end interface ma97_lmultiply

   interface ma97_sparse_fwd_solve
      module procedure ma97_sparse_fwd_solve_double
   end interface ma97_sparse_fwd_solve

   interface ma97_get_n__
      module procedure ma97_get_n_double
   end interface ma97_get_n__

   interface ma97_get_nz__
      module procedure ma97_get_nz_double
   end interface ma97_get_nz__
   
   ! Make smalloc generic
   interface smalloc
      module procedure real_alloc, int_alloc
   end interface smalloc

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!****************************************************************************
!
! Analyse phase.
! Matrix entered in CSC format (lower triangle).
! The user optionally inputs the pivot order. If not, HSL_MC68 called. 
! Structure is then expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
! There is no checking of the user's data if check = .false.
! Otherwise, HSL_MC69 used to clean data.
!
subroutine analyse_double(check, n, ptr, row, akeep, control, info, &
      order, val)
   logical, intent(in) :: check ! if set to true, matrix data is checked
     ! and cleaned data stored in akeep
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: row(:) ! row indices of lower triangular part
   integer, intent(in) :: ptr(:) ! col pointers for lower triangular part
   type(ma97_akeep), intent(out) :: akeep ! See derived-type declaration
   type(ma97_control), intent(in) :: control ! See derived-type declaration
   type(ma97_info), intent(out) :: info      ! See derived-type declaration
   integer, optional, intent(inout) :: order(:)
     ! Must be present and set on entry if control%ordering = 0.
     ! If control%ordering = 0 and i is used to index a variable, |order(i)|
     ! must hold its position in the pivot sequence. If a 1x1 pivot i is
     ! required, the user must set order(i)>0. If a 2x2 pivot involving
     ! variables i and j is required, the user must set
     ! order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
     ! If i is not used to index a variable,
     ! order(i) must be set to zero.
     ! On exit, holds the pivot order to be used by factorization.
     ! Note: this input is consistent with our out-of-core solvers.
     !!!!! Note: although we allow 2x2 pivots to be input, we actually ignore 
     ! the signs (we reset signs of order after call to hsl_mc68 or hsl_mc80)
   real(wp), optional, intent(in) :: val(:) ! must be present
     ! if a matching-based elimination ordering is required
     ! (control%ordering 7 or 8).
     ! If present,  val(k) must hold the value of the entry in row(k).

   character(50)  :: context      ! Procedure name (used when printing).
   integer :: flag69       ! error flag for hsl_mc69
   integer :: mp           ! stream number for diagnostic messages
   integer :: i,j
   integer :: nout         ! stream for errors
   integer :: nout1        ! stream for warnings
   integer :: nz           ! entries in expanded matrix
   integer :: ord80        ! controls ordering of compressed matrix within mc80
   integer :: st           ! stat parameter

   integer, dimension(:), allocatable :: order2
   integer, dimension(:), allocatable :: perm
   integer, dimension(:), allocatable :: ptr2 ! col. pointers for expanded mat
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix

   ! the following are only used for matching-based orderings
   real(wp), dimension(:), allocatable :: val_clean ! cleaned values if
     ! val is present and checking is required. 
   real(wp), dimension(:), allocatable :: val2 ! expanded matrix if
     ! val is present.

   type(mc80_control) :: control80
   type(mc80_info) :: info80

   ! Initialise
   context = 'ma97_analyse'
   call ma97_free(akeep)
   info%flag = 0
   info%matrix_missing_diag = 0
   info%matrix_outrange = 0
   info%matrix_dup = 0
   info%matrix_rank = n
   info%maxdepth = 0
   info%num_sup = 0
   info%ordering = 0
   info%stat = 0

   ! Set stream numbers
   mp = control%unit_diagnostics
   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   nout1 = control%unit_warning
   if (control%print_level < 0) nout1 = -1

   ! Print status on entry
   if (control%print_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to ma97_analyse:'
     write (mp,'(a,i15)') ' control%print_level       =  ', &
        control%print_level
     write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
     write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
     write (mp,'(a,i15)') ' control%nemin             =  ',control%nemin
     write (mp,'(a,i15)') ' control%ordering          =  ',control%ordering
     write (mp,'(a,i15)') ' n                         =  ',n
   end if

   akeep%check = check
   akeep%n = n
   akeep%flag = 0

   ! Checking of matrix data
   if (n < 0) then
      info%flag = MA97_ERROR_A_N_OOR
      call ma97_print_flag(context,nout,info%flag)
      akeep%flag = info%flag
      return
   end if

   if (n .eq. 0) then
      akeep%nnodes = 0
      allocate(akeep%sptr(0), stat=st) ! used to check if analyse has been run
      if (st .ne. 0) go to 490
      akeep%matrix_dup = 0
      akeep%matrix_missing_diag = 0
      akeep%matrix_outrange = 0
      akeep%maxdepth = 0
      akeep%num_sup = 0
      akeep%ordering = 0
      return
   end if

   ! check control%ordering has a valid value
   if (control%ordering < 0 .or. control%ordering > 8) then
      info%flag = MA97_ERROR_ORDER
      call ma97_print_flag(context,nout,info%flag)
      akeep%flag = info%flag
      return
   end if

   ! check val present when expected
   if (control%ordering.eq.7 .or. control%ordering.eq.8) then
     if (.not.present(val)) then
        info%flag = MA97_ERROR_VAL
        call ma97_print_flag(context,nout,info%flag)
        akeep%flag = info%flag
        return
     end if
   end if

   akeep%ne = ptr(n+1)-1

   st = 0
   if (check) then
      deallocate (akeep%ptr,stat=st)
      allocate (akeep%ptr(n+1),stat=st)
      if (st .ne. 0) go to 490

      if (present(val)) then
         call mc69_cscl_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
            akeep%ptr, akeep%row, flag69, val_in=val, val_out=val_clean,   &
            lmap=akeep%lmap, map=akeep%map,  &
            noor=info%matrix_outrange, ndup=info%matrix_dup)
      else
         call mc69_cscl_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, &
            akeep%ptr, akeep%row, flag69, lmap=akeep%lmap, map=akeep%map,  &
            noor=info%matrix_outrange, ndup=info%matrix_dup)
      end if
      ! Check for errors
      if (flag69 < 0) then
         if (flag69 .eq. -1) info%flag  = MA97_ERROR_ALLOCATION
         if (flag69 .eq. -5) info%flag  = MA97_ERROR_A_PTR
         if (flag69 .eq. -6) info%flag  = MA97_ERROR_A_PTR
         if (flag69 .eq. -10) info%flag = MA97_ERROR_A_ALL_OOR
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      end if

      ! Check whether warning needs to be raised
      ! Note: same numbering of positive flags as in mc69
      if (flag69 > 0) then
         info%flag = flag69
         call ma97_print_flag(context,nout1,info%flag)
      end if
      nz = akeep%ptr(n+1) - 1
   else
      nz = ptr(n+1)-1
   end if

   ! If the pivot order is not supplied, we need to compute an order.
   ! Otherwise, we check the supplied order.

   deallocate (akeep%invp,stat=st)
   allocate (akeep%invp(n),perm(n),order2(n),ptr2(n+1),row2(2*nz),stat=st)
   if (st .ne. 0) go to 490
   if(control%ordering.eq.7 .or. control%ordering.eq.8) then
      allocate(akeep%scaling(n), val2(2*nz), stat=st)
      if (st .ne. 0) go to 490
   end if

   select case(control%ordering)
   case(0)
      if (.not.present(order)) then
         ! we have an error since user should have supplied the order
         info%flag = MA97_ERROR_ORDER
         akeep%flag = info%flag
         call ma97_print_flag(context,nout,info%flag)
         return
      end if
      call check_order(n,order,akeep%invp,perm,akeep,control,info)
      if (info%flag < 0) go to 490
      order2(1:n) = order(1:n)
      if (check) then
         call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
      else
         call expand_pattern(n, nz, ptr, row, ptr2, row2)
      end if
      info%ordering = 0
   case(7:8)
      ! matching-based ordering required
      ! Expand the matrix as more efficient to do it and then
      ! call hsl_mc80 with full matrix supplied

      if (check) then
         call expand_matrix(n, nz, akeep%ptr, akeep%row, val_clean, ptr2, &
            row2, val2)
         deallocate (val_clean,stat=st)
      else
         call expand_matrix(n, nz, ptr, row, val, ptr2, row2, val2)
      end if

      if (control%ordering .eq. 7) ord80 = 1
      if (control%ordering .eq. 8) ord80 = 3

      !control80%unmatched_last = .true.
      !control80%unmatched_scale_zero = .true.

      call mc80_order_full(ord80, n, ptr2, row2, val2, order2, control80, &
         info80, scale=akeep%scaling)

      select case(info80%flag)
      case(0)
         ! Success; do nothing
      case(1)
         ! singularity warning required
         info%flag = MA97_WARNING_ANAL_SINGULAR
         call ma97_print_flag(context,nout1,info%flag)
      case(-1)
         info%flag = MA97_ERROR_ALLOCATION
         info%stat = info80%stat
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case(-4)
         info%flag = MA97_ERROR_MC68
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case(-6)
         info%flag = MA97_ERROR_NO_METIS
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case default
         info%flag = MA97_ERROR_UNKNOWN
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      end select

      info%matrix_rank = info80%compress_rank
      info%ordering = control%ordering
      ! we do not allow input of 2x2 pivots and so we get rid of negative signs
      do i = 1,n
         j = abs(order2(i))
         perm(i) = j
      end do
      deallocate (val2,stat=st)
   case default
      ! ordering computed using hsl_mc68 and then pattern of matrix expanded
      if (check) then
         call compute_order(n,nz,akeep%ptr,akeep%row,order2,akeep%invp,perm, &
            control,info)
         call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
      else
         call compute_order(n,nz,ptr,row,order2,akeep%invp,perm, &
            control,info)
         call expand_pattern(n, nz, ptr, row, ptr2, row2)
      end if
      if (info%flag < 0) go to 490
   end select


   ! perform rest of analyse
   if (check) then
      call analyse_phase(n, akeep%ptr, akeep%row, ptr2, row2, order2,  &
         akeep%invp, perm, akeep, control, info)
   else
      call analyse_phase(n, ptr, row, ptr2, row2, order2, akeep%invp, perm, &
         akeep, control, info)
   end if

   if (present(order)) order(1:n) = abs(order2(1:n))

   490 continue
   info%stat = st
   if (info%stat .ne. 0) then
      info%flag = MA97_ERROR_ALLOCATION
      call ma97_print_flag(context,nout,info%flag,st=info%stat)
   end if
   akeep%flag = info%flag

end subroutine analyse_double

!****************************************************************************
!
! Analyse phase.
! Matrix entered in coordinate format.
! HSL_MC69 is used to convert the data to CSC format.
! The user optionally inputs the pivot order. If not, HSL_MC68 called. 
! Structure is then expanded.
! Supervariables are computed
! and then the assembly tree is constructed and the data structures
! required by the factorization are set up.
!
subroutine ma97_analyse_coord_double(n, ne, row, col, akeep, control, &
      info, order, val)
   integer, intent(in) :: n ! order of A
   integer, intent(in) :: ne ! entries to be input by user
   integer, intent(in) :: row(:) ! row indices
   integer, intent(in) :: col(:) ! col indices
   type(ma97_akeep), intent(out) :: akeep ! See derived-type declaration
   type(ma97_control), intent(in) :: control ! See derived-type declaration
   type(ma97_info), intent(out) :: info      ! See derived-type declaration
   integer, intent(inout), optional  :: order(:)
      ! Must be present and set on entry if control%ordering = 0 
      ! i is used to index a variable, |order(i)| must
      ! hold its position in the pivot sequence. If a 1x1 pivot i is required,
      ! the user must set order(i)>0. If a 2x2 pivot involving variables
      !  i and j is required, the user must set
      !  order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
      !  If i is not used to index a variable,
      !  order(i) must be set to zero.
      !  On exit, holds the pivot order to be used by factorization.
      !!!!! Note: although we allow 2x2 pivots to be input, we actually ignore 
      ! the signs (we reset signs of order after call to hsl_mc68 or hsl_mc80)
   real(wp), optional, intent(in) :: val(:) ! must be present
     ! if a matching-based elimination ordering is required 
     ! (control%ordering 7 or 8).
     ! If present, val(k) must hold value of entry in row(k) and col(k).

    integer, dimension(:), allocatable :: perm
      ! Allocated to have size n.

   integer, dimension(:), allocatable :: ptr2 ! col. pointers for expanded mat
   integer, dimension(:), allocatable :: row2 ! row indices for expanded matrix
   integer, dimension(:), allocatable :: order2 ! pivot order

   type(mc80_control) :: control80
   type(mc80_info) :: info80

   real(wp), dimension(:), allocatable :: val_clean ! cleaned values if
     ! val is present.
   real(wp), dimension(:), allocatable :: val2 !expanded matrix (val present)

   character(50)  :: context      ! Procedure name (used when printing).
   integer :: flag69       ! error flag for hsl_mc69
   integer :: i
   integer :: mp           ! stream number for diagnostic messages
   integer :: nout         ! stream for errors
   integer :: nout1        ! stream for warnings
   integer :: nz           ! entries in expanded matrix
   integer :: ord80        ! controls ordering of compressed matrix within mc80
   integer :: st           ! stat parameter

   ! Initialise
   context = 'ma97_analyse_coord'
   call ma97_free(akeep)
   info%flag = 0
   info%matrix_missing_diag = 0
   info%matrix_outrange = 0
   info%matrix_dup = 0
   info%matrix_rank = n
   info%maxdepth = 0
   info%num_sup = 0
   info%ordering = 0
   info%stat = 0

   ! Set stream numbers
   mp = control%unit_diagnostics
   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   nout1 = control%unit_warning
   if (control%print_level < 0) nout1 = -1

   ! Output status on entry
   if (control%print_level>=1 .and. mp>=0) then
     write (mp,'(/a)') ' On entry to ma97_analyse_coord:'
     write (mp,'(a,i15)') ' control%print_level       =  ', &
        control%print_level
     write (mp,'(a,i15)') ' control%unit_diagnostics  =  ',mp
     write (mp,'(a,i15)') ' control%unit_error        =  ',control%unit_error
     write (mp,'(a,i15)') ' control%unit_warning      =  ',control%unit_warning
     write (mp,'(a,i15)') ' control%nemin             =  ',control%nemin
     write (mp,'(a,i15)') ' control%ordering          =  ',control%ordering
     write (mp,'(a,i15)') ' n                         =  ',n
     write (mp,'(a,i15)') ' ne                        =  ',ne
   end if

   akeep%check = .true.
   akeep%n = n
   akeep%ne = ne
   akeep%flag = 0

   !
   ! Checking of matrix data
   !
   if (n < 0 .or. ne < 0) then
      info%flag = MA97_ERROR_A_N_OOR
      akeep%flag = info%flag
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (n .eq. 0) then
      akeep%nnodes = 0
      allocate(akeep%sptr(0), stat=st) ! used to check if analyse has been run
      if (st .ne. 0) go to 490
      akeep%matrix_dup = 0
      akeep%matrix_missing_diag = 0
      akeep%matrix_outrange = 0
      akeep%maxdepth = 0
      akeep%num_sup = 0
      akeep%ordering = 0
      return
   end if

   ! check control%ordering has a valid value
   if (control%ordering < 0 .or. control%ordering > 8) then
      info%flag = MA97_ERROR_ORDER
      call ma97_print_flag(context,nout,info%flag)
      akeep%flag = info%flag
      return
   end if

   ! check val present when expected
   if (control%ordering.eq.7 .or. control%ordering.eq.8) then
     if (.not.present(val)) then
        info%flag = MA97_ERROR_VAL
        call ma97_print_flag(context,nout,info%flag)
        akeep%flag = info%flag
        return
     end if
   end if

   st = 0
   deallocate (akeep%ptr,stat=st)
   allocate (akeep%ptr(n+1),stat=st)
   if (st .ne. 0) go to 490

   if (present(val)) then
      call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, row, col, &
         akeep%ptr, akeep%row, flag69, val_in=val, val_out=val_clean,        &
         lmap=akeep%lmap, map=akeep%map, &
         noor=info%matrix_outrange,  ndup=info%matrix_dup)
   else
      call mc69_coord_convert(HSL_MATRIX_REAL_SYM_INDEF, n, n, ne, row, col, &
         akeep%ptr, akeep%row, flag69, lmap=akeep%lmap, map=akeep%map, &
         noor=info%matrix_outrange,  ndup=info%matrix_dup)
   end if

   ! Check for errors
   if (flag69 < 0) then
      if (flag69 .eq. -1)  info%flag = MA97_ERROR_ALLOCATION
      if (flag69 .eq. -10) info%flag = MA97_ERROR_A_ALL_OOR
      call ma97_print_flag(context,nout,info%flag)
      akeep%flag = info%flag
      return
   end if

   ! Check whether warning needs to be raised
   ! Note: same numbering of positive flags as in mc69
   if (flag69 > 0) then
      info%flag = flag69
      call ma97_print_flag(context,nout1,info%flag)
      akeep%flag = info%flag
   end if

   nz = akeep%ptr(n+1) - 1

   ! If the pivot order is not supplied, we need to compute an order
   ! here, before we expand the matrix structure.
   ! Otherwise, we must check the supplied order.

   deallocate(akeep%invp, stat=st)
   allocate (akeep%invp(n),perm(n),order2(n),ptr2(n+1),row2(2*nz),stat=st)
   if (st .ne. 0) go to 490
   if(control%ordering.eq.7 .or. control%ordering.eq.8) then
      allocate(val2(2*nz),akeep%scaling(n),stat=st)
      if (st .ne. 0) go to 490
   end if

   select case(control%ordering)
   case(0)
      if (.not.present(order)) then
         ! we have an error since user should have supplied the order
         info%flag = MA97_ERROR_ORDER
         akeep%flag = info%flag
         call ma97_print_flag(context,nout,info%flag)
         return
      end if
      call check_order(n,order,akeep%invp,perm,akeep,control,info)
      if (info%flag < 0) go to 490
      order2(1:n) = order(1:n)
      call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
      info%ordering = 0

   case(7:8)
      ! matching-based ordering required

      call expand_matrix(n, nz, akeep%ptr, akeep%row, val_clean, ptr2, row2, &
         val2)
      deallocate (val_clean,stat=st)

      if (control%ordering .eq. 7) ord80 = 1
      if (control%ordering .eq. 8) ord80 = 3

      call mc80_order_full(ord80,n,ptr2,row2,val2,order2,control80,info80,&
         scale=akeep%scaling)

      select case(info80%flag)
      case(0)
         ! Success; do nothing
      case(1)
         ! singularity warning required
         info%flag = MA97_WARNING_ANAL_SINGULAR
         call ma97_print_flag(context,nout1,info%flag)
      case(-1)
         info%flag = MA97_ERROR_ALLOCATION
         info%stat = info80%stat
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case(-4)
         info%flag = MA97_ERROR_MC68
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case(-6)
         info%flag = MA97_ERROR_NO_METIS
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      case default
         info%flag = MA97_ERROR_UNKNOWN
         call ma97_print_flag(context,nout,info%flag)
         akeep%flag = info%flag
         return
      end select

      info%matrix_rank = info80%compress_rank
      info%ordering = control%ordering
      do i = 1,n
         perm(i) = abs(order2(i))
      end do
      deallocate (val2,stat=st)

   case default
      call compute_order(n,nz,akeep%ptr,akeep%row,order2,akeep%invp,perm, &
         control,info)
      if (info%flag < 0) go to 490
      call expand_pattern(n, nz, akeep%ptr, akeep%row, ptr2, row2)
   end select


   ! we now have the expanded structure held using ptr2, row2
   ! and we are ready to get on with the analyse phase.
   call analyse_phase(n, akeep%ptr, akeep%row, ptr2, row2, order2,  &
      akeep%invp, perm, akeep, control, info)
   if (info%flag < 0) go to 490

   if (present(order)) order(1:n) = abs(order2(1:n))

   490 continue
   info%stat = st
   if (info%stat .ne. 0) then
      info%flag = MA97_ERROR_ALLOCATION
      call ma97_print_flag(context,nout,info%flag,st=info%stat)
   end if
   akeep%flag = info%flag
    
end subroutine ma97_analyse_coord_double

!****************************************************************************

!
! Given lower triangular part of A held in row and ptr, expand to
! upper and lower triangular parts (pattern only). No checks.
!
! Note: we do not use hsl_mc34 here to expand A since, if we did, we would
! an extra copy of the lower triangle into the full structure before
! calling hsl_mc34
!
subroutine expand_pattern(n,nz,ptr,row,aptr,arow)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: nz
   integer, intent(in) :: ptr(n+1)
   integer, intent(in) :: row(nz)
   integer, intent(out) :: aptr(n+1)
   integer, intent(out) :: arow(2*nz)

   integer :: i,j,k

   ! Set aptr(j) to hold no. nonzeros in column j
   aptr(:) = 0
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         aptr(i) = aptr(i) + 1
         if (j.eq.i) cycle
         aptr(j) = aptr(j) + 1
      end do
   end do

   ! Set aptr(j) to point to where row indices will end in arow
   do j = 2, n
      aptr(j) = aptr(j-1) + aptr(j)
   end do
   aptr(n+1) = aptr(n) + 1

   ! Fill arow and aptr
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         arow(aptr(i)) = j
         aptr(i) = aptr(i) - 1
         if (j.eq.i) cycle
         arow(aptr(j)) = i
         aptr(j) = aptr(j) - 1
      end do
   end do
   do j = 1,n
      aptr(j) = aptr(j) + 1
   end do
end subroutine expand_pattern

!****************************************************************************
!
! Given lower triangular part of A held in row, val and ptr, expand to
! upper and lower triangular parts.

subroutine expand_matrix(n,nz,ptr,row,val,aptr,arow,aval)

   integer, intent(in)   :: n ! order of system
   integer, intent(in)   :: nz
   integer, intent(in)   :: ptr(n+1)
   integer, intent(in)   :: row(nz)
   real(wp), intent(in)  :: val(nz)
   integer, intent(out)  :: aptr(n+1)
   integer, intent(out)  :: arow(2*nz)
   real(wp), intent(out) :: aval(2*nz)

   integer :: i,j,k,ipos,jpos
   real(wp) :: atemp

   ! Set aptr(j) to hold no. nonzeros in column j
   aptr(:) = 0
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         aptr(i) = aptr(i) + 1
         if (j.eq.i) cycle
         aptr(j) = aptr(j) + 1
      end do
   end do

   ! Set aptr(j) to point to where row indices will end in arow
   do j = 2, n
      aptr(j) = aptr(j-1) + aptr(j)
   end do
   aptr(n+1) = aptr(n) + 1

   ! Fill arow, aval and aptr
   do j = 1, n
      do k = ptr(j), ptr(j+1) - 1
         i = row(k)
         atemp = val(k)
         ipos = aptr(i)
         arow(ipos) = j
         aval(ipos) = atemp
         aptr(i) = ipos - 1
         if (j.eq.i) cycle
         jpos = aptr(j)
         arow(jpos) = i
         aval(jpos) = atemp
         aptr(j) = jpos - 1
      end do
   end do
   do j = 1,n
      aptr(j) = aptr(j) + 1
   end do

end subroutine expand_matrix


!****************************************************************************
!
! This routine requires the LOWER triangular part of A
! to be held in CSC format.
! The user has supplied a pivot order and this routine checks it is OK
! and returns an error if not. Also sets perm, invp.
!
subroutine check_order(n, order, invp, perm, akeep, control, info)
    integer, intent(in) :: n ! order of system
    integer, intent(inout) :: order(:)
      ! If i is used to index a variable, |order(i)| must
      ! hold its position in the pivot sequence. If 1x1 pivot i required,
      ! the user must set order(i)>0. If a 2x2 pivot involving variables
      ! i and j is required, the user must set
      ! order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
      ! If i is not used to index a variable, order(i) must be set to zero.
      ! !!!! In this version, signs are reset to positive value
   integer, intent(out) :: invp(n)
      ! Used to check order and then holds inverse of perm.
   integer, intent(out) :: perm(n)
   type (ma97_akeep), intent(inout) :: akeep
   type (ma97_control), intent(in) :: control
   type (ma97_info), intent(inout) :: info

   character(50)  :: context ! Procedure name (used when printing).

   integer :: i, j
   integer :: k, l
   integer :: nout  ! stream for error messages
   integer :: num_null ! number of unused variables

   context = 'ma97_analyse'
   nout = control%unit_error
   if (control%print_level < 0) nout = -1

   if (size(order) < n) then
      ! Order is too short
      info%flag = MA97_ERROR_ORDER
      akeep%flag = info%flag
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   ! initialise
   num_null = 0
   invp(:) = 0

   do i = 1,n
      order(i) = abs(order(i))
   end do
     
   ! Check user-supplied order and copy the absolute values to invp.
   ! Also add up number of variables that are not used (null rows)
   do i = 1, n
      j = order(i)
      if (j .eq. 0) then
         num_null = num_null + 1
         cycle
      end if
      if (j > n) exit
      if (invp(j) .ne. 0) exit ! Duplicate found
      invp(j) = i
   end do
   if (i-1 .ne. n) then
      info%flag = MA97_ERROR_ORDER
      akeep%flag = info%flag
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   ! ensure invp is set for null rows/cols (want abs(invp) to hold
   ! a permutation)
   if (num_null > 0) then
      k = 1
      do i = 1,n
         j = order(i)
         if (j.ne.0) cycle
         do l = k,n
            if (invp(l) .eq. 0) exit
         end do
         invp(l) = i
         k = l + 1
      end do
   end if

   ! set perm to be inverse of invp (=abs(order) if num_null = 0)
   do i = 1, n
      j = abs(invp(i))
      perm(j) = i
   end do
 
end subroutine check_order

!****************************************************************************
!
! This routine requires the LOWER triangular part of A
! to be held in CSC format.
! This routines computes a pivot order.
!
subroutine compute_order(n, ne, ptr, row, order, invp, perm,  &
      control, info)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ne ! entries in lower triangle of A
   integer, intent(in) :: row(ne) ! row indices 
   integer, intent(in) :: ptr(n+1) ! col pointers
   integer, intent(out) :: order(:) ! set by HSL_MC68
   integer, intent(out) :: invp(n)
      ! Used to check order and then holds inverse of perm.
   integer, intent(out) :: perm(n)
   type (ma97_control), intent(in) :: control
   type (ma97_info), intent(inout) :: info

   ! Following parameters are used for heuristic choice between AMD and
   ! MeTiS. See report RAL-TR-2006-001:
   ! "Towards an automatic ordering for a symmetric sparse direct solver"
   ! I.S. Duff and J.A. Scott
   integer, parameter :: bigN = int(1e5)
   real(wp), parameter :: c(3) = (/ 3.0, 1.8, 10.0 /)

   type (mc68_control) :: control68

   character(50)  :: context ! Procedure name (used when printing).

   integer :: nout  ! stream for error messages
   integer :: nout1 ! stream for warnings
   integer :: ord68
   integer :: flag, st

   ! Variables used in automatic determination of ordering method
   logical :: oxo
   integer(long) :: amd_nzl, metis_nzl
   integer :: i, j, sz11, sz22
   integer, dimension(:), allocatable :: order2, invp2, perm2

   context = 'ma97_analyse'
   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   nout1 = control%unit_warning
   if (control%print_level < 0) nout1 = -1

   ! switch off hsl_mc68 printing
   control68%lp = -1
   control68%wp = -1
   control68%mp = -1
   control68%print_level = -1

   select case(control%ordering)
   case(1:4) ! As per equivalent mc68 option: AMDD, AMD, MeTiS, MA47
      ord68 = control%ordering
      call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
         control68, flag, info%stat, info%flag68)
   case(5) ! Automatic choice between AMD and MeTiS (parallel)
      ! FIXME: This heuristic could do with further validation
      if(n.gt.bigN) then
         ! n is moderatly large, just use MeTiS if possible
         ord68 = 3 ! MeTiS
         call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
            control68, flag, info%stat, info%flag68)
         if(flag.eq.MA97_ERROR_NO_METIS) then
            ord68 = 1 ! AMD
            call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
               control68, flag, info%stat, info%flag68)
         end if
      else ! n is small, use same heuristic as in serial
         ! Check if we have an OXO matrix:
         ! Count size of a possible (2,2) block
         sz22 = 0
         do i = n, 1, -1
            if(ptr(i).ne.ptr(i+1)) exit
            sz22 = sz22 + 1
         end do
         oxo = (sz22.gt.0)
         sz11 = n - sz22
         ! Check if matching (1,1) block is also empty
         oxo22a: do i = 1, sz11
            do j = ptr(i), ptr(i+1)-1
               if(row(j).lt.sz11) then
                  oxo = .false.
                  exit oxo22a
               end if
            end do
         end do oxo22a
         if(oxo) then
            if(sz11 .gt. c(2)*sz22) then
               ord68 = 1 ! AMD
               call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                  control68, flag, info%stat, info%flag68)
            else
               ord68 = 3 ! MeTiS
               call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                  control68, flag, info%stat, info%flag68)
               if(flag.eq.MA97_ERROR_NO_METIS) then
                  ord68 = 1 ! AMD
                  call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                     control68, flag, info%stat, info%flag68)
               end if
            end if
         else
            ord68 = 1 ! AMD
            call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
               control68, flag, info%stat, info%flag68)
            ! Find predicted number of entries in factors
            amd_nzl = get_nzl(n, ptr, row, perm, st)
            if(st.ne.0) goto 10
            if(amd_nzl .gt. c(3)*ne) then
               ! Handle errors from call for AMD
               if(flag.lt.0) then
                  info%flag = flag
                  if (info%flag .ne. MA97_ERROR_ALLOCATION) &
                     call ma97_print_flag(context,nout,info%flag,st=info%stat)
                  return
               elseif(flag.gt.0) then
                  info%flag = flag
                  call ma97_print_flag(context,nout1,info%flag) 
               end if
               ! Allocate arrays to store MeTiS ordering temporarily
               allocate(order2(n), invp2(n), perm2(n), stat=st)
               if(st.ne.0) goto 10
               ord68 = 3 ! MeTiS
               call mc68_wrap(ord68, n, ne, ptr, row, order2, invp2, perm2,  &
                  control68, flag, info%stat, info%flag68)
               if(flag.eq.MA97_ERROR_NO_METIS) then
                  ord68 = 1 ! AMD
                  flag = 0; info%flag68 = 0; info%stat = 0
               elseif(flag.ge.0) then
                  ! Find predicted number of entries in factors
                  metis_nzl = get_nzl(n, ptr, row, perm2, st)
                  if(st.ne.0) goto 10
                  if(metis_nzl .lt. amd_nzl) then
                     ! Use MeTiS ordering
                     ord68 = 3 ! MeTiS
                     order(:) = order2(:)
                     invp(:) = invp2(:)
                     perm(:) = perm2(:)
                  else
                     ! Use AMD ordering
                     ord68 = 1 ! AMD
                  end if
               end if
            end if
         end if
      end if
   case(6) ! Automatic choice between AMD and MeTiS (serial)
      ! For more details, see RAL-TR-2006-001:
      ! "Towards an automatic ordering for a symmetric sparse direct solver"
      ! I.S. Duff and J.A. Scott
      if(n.gt.bigN) then
         if(ne .lt. c(1)*n) then
            ord68 = 1 ! AMD
            call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
               control68, flag, info%stat, info%flag68)
         else
            ord68 = 3 ! MeTiS
            call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
               control68, flag, info%stat, info%flag68)
            if(flag.eq.MA97_ERROR_NO_METIS) then
               ord68 = 1 ! AMD
               call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                  control68, flag, info%stat, info%flag68)
            end if
         end if
      else
         ! Check if we have an OXO matrix:
         ! Count size of a possible (2,2) block
         sz22 = 0
         do i = n, 1, -1
            if(ptr(i).ne.ptr(i+1)) exit
            sz22 = sz22 + 1
         end do
         oxo = (sz22.gt.0)
         sz11 = n - sz22
         ! Check if matching (1,1) block is also empty
         oxo22b: do i = 1, sz11
            do j = ptr(i), ptr(i+1)-1
               if(row(j).lt.sz11) then
                  oxo = .false.
                  exit oxo22b
               end if
            end do
         end do oxo22b
         if(oxo) then
            if(sz11 .gt. c(2)*sz22) then
               ord68 = 1 ! AMD
               call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                  control68, flag, info%stat, info%flag68)
            else
               ord68 = 3 ! MeTiS
               call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                  control68, flag, info%stat, info%flag68)
               if(flag.eq.MA97_ERROR_NO_METIS) then
                  ord68 = 1 ! AMD
                  call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
                     control68, flag, info%stat, info%flag68)
               end if
            end if
         else
            ord68 = 1 ! AMD
            call mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
               control68, flag, info%stat, info%flag68)
            ! Find predicted number of entries in factors
            amd_nzl = get_nzl(n, ptr, row, perm, st)
            if(st.ne.0) goto 10
            if(amd_nzl .gt. c(3)*ne) then
               ! Handle errors from call for AMD
               if(flag.lt.0) then
                  info%flag = flag
                  if (info%flag .ne. MA97_ERROR_ALLOCATION) &
                     call ma97_print_flag(context,nout,info%flag,st=info%stat)
                  return
               elseif(flag.gt.0) then
                  info%flag = flag
                  call ma97_print_flag(context,nout1,info%flag) 
               end if
               ! Allocate arrays to store MeTiS ordering temporarily
               allocate(order2(n), invp2(n), perm2(n), stat=st)
               if(st.ne.0) goto 10
               ord68 = 3 ! MeTiS
               call mc68_wrap(ord68, n, ne, ptr, row, order2, invp2, perm2,  &
                  control68, flag, info%stat, info%flag68)
               if(flag.eq.MA97_ERROR_NO_METIS) then
                  ord68 = 1 ! AMD
                  flag = 0; info%flag68 = 0; info%stat = 0
               elseif(flag.ge.0) then
                  ! Find predicted number of entries in factors
                  metis_nzl = get_nzl(n, ptr, row, perm2, st)
                  if(st.ne.0) goto 10
                  if(metis_nzl .lt. amd_nzl) then
                     ! Use MeTiS ordering
                     ord68 = 3 ! MeTiS
                     order(:) = order2(:)
                     invp(:) = invp2(:)
                     perm(:) = perm2(:)
                  else
                     ! Use AMD ordering
                     ord68 = 1 ! AMD
                  end if
               end if
            end if
         end if
      end if
   case default
      info%flag = MA97_ERROR_ORDER
      call ma97_print_flag(context,nout,info%flag,st=info%stat)
      return
   end select
   info%ordering = ord68

   if(flag.lt.0) then
      info%flag = flag
      if (info%flag .ne. MA97_ERROR_ALLOCATION) &
         call ma97_print_flag(context,nout,info%flag,st=info%stat)
      return
   elseif(flag.gt.0) then
      info%flag = flag
      call ma97_print_flag(context,nout1,info%flag) 
   end if

   return

   10 continue
   info%flag = MA97_ERROR_ALLOCATION
   info%stat = st
   call ma97_print_flag(context,nout,info%flag,st=info%stat)
   return
end subroutine compute_order

!****************************************************************************

!
! Following function is based on mc78_analyse, but only goes as far as
! required to calculate the number of entries in L (without amalgamation)
!
integer(long) function get_nzl(n, ptr, row, perm, st)
   integer, intent(in) :: n ! Dimension of system
   integer, dimension(n+1), intent(in) :: ptr ! Column pointers
   integer, dimension(ptr(n+1)-1), intent(in) :: row ! Row indices
   integer, dimension(n), intent(in) :: perm
      ! perm(i) must hold position of i in the pivot sequence. 
   integer, intent(out) :: st

   integer :: realn
   integer :: i
   integer, dimension(:), allocatable :: cc, parent, perm2, invp
   integer, dimension(:), allocatable :: ptr2, row2

   get_nzl = huge(0)

   ! Allocate memory
   allocate(invp(n+1), perm2(n), parent(n), cc(n+1), ptr2(n+1), &
      row2(2*ptr(n+1)), stat=st)
   if(st.ne.0) return

   ! Expand matrix
   ptr2(1:n+1) = ptr(1:n+1)
   row2(1:ptr(n+1)-1) = row(1:ptr(n+1)-1)
   call mc34_expand(n, row2, ptr2, invp) ! invp used as workspace

   ! Initialise inverse permutation
   perm2(1:n) = perm(1:n)
   do i = 1, n
      invp(perm(i)) = i
   end do

   ! Build elimination tree
   call mc78_etree(n, ptr2, row2, perm2, invp, parent, st)
   if(st.ne.0) return

   ! Postorder tree (modifies perm!)
   call mc78_postorder(n, realn, ptr2, perm2, invp, parent, st)
   if(st.ne.0) return

   ! Determine column counts
   call mc78_col_counts(n, ptr2, row2, perm2, invp, parent, cc, st)
   if(st.ne.0) return

   ! Calculate stats
   do i = 1, n+1
      invp(i) = i
   end do
   call mc78_stats(n, invp, cc, nfact=get_nzl)
end function get_nzl

!****************************************************************************

! This subroutine wraps a call to mc68_order with the given ord68 ordering
! to minimise error checking and other post-processing repetition.
subroutine mc68_wrap(ord68, n, ne, ptr, row, order, invp, perm,  &
      control68, flag, stat, flag68)
   integer, intent(in) :: ord68
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ne ! entries in lower triangle of A
   integer, intent(in) :: row(ne) ! row indices 
   integer, intent(in) :: ptr(n+1) ! col pointers
   integer, intent(out) :: order(:) ! set by HSL_MC68
   integer, intent(out) :: invp(n)
      ! Used to check order and then holds inverse of perm.
   integer, intent(out) :: perm(n)
   type (mc68_control), intent(in) :: control68
   integer, intent(out) :: flag
   integer, intent(out) :: stat
   integer, intent(out) :: flag68

   ! mc68 derived types
   type (MC68_info) :: info68

   integer :: i, j
   integer :: k
   integer :: l
   integer :: num_null ! number of unused variables

   flag = 0
   flag68 = 0
   stat = 0

   call mc68_order(ord68,n,ptr,row,order,control68,info68)

   if (info68%flag < 0) then
      select case(info68%flag)
      case(-1)
         flag = MA97_ERROR_ALLOCATION
         stat = info68%stat
      case(-4)
         flag = MA97_ERROR_ORDER
      case(-5)
         flag = MA97_ERROR_NO_METIS
      case default
         flag = MA97_ERROR_MC68
         flag68 = info68%flag
      end select
      return
   elseif (info68%flag > 1) then
      flag = MA97_WARNING_ANAL_SINGULAR
   end if

   if (ord68 .eq. 4) then
      !!! the pivot sequence will contain 2x2 pivots. In this
      !   version we are going to ignore 2x2 pivots 
      do i = 1,n
         order(i) = abs(order(i))
      end do
   end if

   ! Add up number of variables that are not used (null rows)

   !!!! HSL_MC68 spec sheet says that order(i) = 0 if row i is null
   ! but is this correct? Looks as though order contains a permutation
   ! (at least if min deg-type algorithm used) so maybe don't
   ! need the checks on order(i) being 0??

   num_null = 0
   invp(:) = 0

   do i = 1, n
      j = order(i)
      if (j .eq. 0) then
         ! Variable i is not used
         num_null = num_null + 1
         cycle
      end if
      invp(j) = i
   end do

   ! ensure invp is set for null rows/cols (want abs(invp) to hold
   ! a permutation)
   if (num_null > 0) then
      k = 1
      do i = 1,n
         j = order(i)
         if (j.ne.0) cycle
         do l = k,n
            if (invp(l) .eq. 0) exit
         end do
         invp(l) = i
         k = l + 1
      end do
   end if

   ! set perm to be inverse of invp (=abs(order) if num_null = 0)
   do i = 1, n
      j = abs(invp(i))
      perm(j) = i
   end do
end subroutine mc68_wrap

!****************************************************************************

!
! This routine requires the LOWER and UPPER triangular parts of A
! to be held in CSC format using ptr2 and row2
! AND lower triangular part held using ptr and row.
!
! On exit from this routine, order is set to order
! input to factorization.
!
subroutine analyse_phase(n, ptr, row, ptr2, row2, order, invp, perm,  &
      akeep, control, info)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! col pointers (lower triangle) 
   integer, intent(in) :: row(ptr(n+1)-1) ! row indices (lower triangle)
   integer, intent(in) :: ptr2(n+1) ! col pointers (whole matrix)
   integer, intent(in) :: row2(ptr2(n+1)-1) ! row indices (whole matrix)
   integer, dimension(n), intent(inout) :: order
      !  On exit, holds the pivot order to be used by factorization.
   integer, dimension(n), intent(out) :: invp 
      ! Work array. Used to hold inverse of order but
      ! is NOT set to inverse for the final order that is returned.
   integer, dimension(n), intent(inout) :: perm
      ! On entry to mc78_analyse, perm(i) holds position of i in pivot 
      ! sequence. On exit from mc78_analyse holds the pivot order 
      ! to be used by ma97_factor. perm is copied to order
      ! before return to user.
   type (ma97_akeep), intent(inout) :: akeep
   type (ma97_control), intent(in) :: control
   type (ma97_info), intent(inout) :: info

   type(mc78_control) :: control78

   character(50)  :: context ! Procedure name (used when printing).
   integer, dimension(:), allocatable :: child_next, child_head ! linked
      ! list for children, used to build akeep%child_ptr, akeep%child_list

   integer  :: blkm, blkn
   integer  :: info78 ! error flag for hsl_mc78
   integer  :: i, j , k
   integer  :: nout, nout1 ! streams for errors and warnings
   integer  :: nz ! ptr(n+1)-1
   integer  :: st

   context = 'ma97_analyse'
   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   nout1 = control%unit_warning
   if (control%print_level < 0) nout1 = -1

   ! Check nemin (a node is merged with its parent if both involve
   ! fewer than nemin eliminations). If out of range, use the default
   control78%nemin = control%nemin
   if (control78%nemin < 1) control78%nemin = nemin_default

   ! switch off printing from hsl_mc78
   control78%unit_error = -1
   control78%unit_warning = -1

   ! we want null/rows columns to raise a warning
   control78%ssa_abort = .false.

   ! switch on sorting within each supernode row list.
   control78%sort = .true.

   ! Uncomment next line to do MA57 style s/n amalgamation
   !control78%heuristic = 3; control78%nemin = 16

   call mc78_analyse(n, ptr2, row2, perm, akeep%nnodes, akeep%sptr, &
      akeep%sparent, akeep%rptr,akeep%rlist,                        &
      control78, info78, stat=st, nfact=akeep%nfactor,              &
      nflops=info%num_flops)
   info%num_factor = akeep%nfactor

   ! only posible non-zero values for error flag are +1 and -1.
   select case(info78)
   case(1)
      if (info%flag .ne. MA97_WARNING_ANAL_SINGULAR) then
         info%flag = MA97_WARNING_ANAL_SINGULAR
         call ma97_print_flag(context,nout1,info%flag)
      end if
   case(-1)
      ! allocation error
      go to 490
   end select

   ! copy perm into order and set invp to hold inverse of order
   do i = 1,n
      j = perm(i)
      order(i) = j
      invp(j) = i
   end do
   ! any unused variables are at the end and so can set order for them
   do j = akeep%sptr(akeep%nnodes+1), n
      i = invp(j)
      order(i) = 0
   end do

   ! Build map from A to L in nptr, nlist
   nz = ptr(n+1) - 1
   allocate(akeep%nptr(n+1), akeep%nlist(2,nz), stat=st)
   if (st .ne. 0) go to 490

   call build_map(n, ptr, row, order, invp, akeep%nnodes, akeep%sptr, &
      akeep%rptr, akeep%rlist, akeep%nptr, akeep%nlist, st)
   if (st .ne. 0) go to 490

   ! Find maxmn and setup levels
   allocate(akeep%level(akeep%nnodes+1), stat=st)
   if (st .ne. 0) go to 490

   akeep%maxmn = 0
   akeep%level(akeep%nnodes+1) = 0
   info%maxfront = 0
   info%maxdepth = 0
   do i = akeep%nnodes, 1, -1
      blkn = akeep%sptr(i+1) - akeep%sptr(i) 
      blkm = int(akeep%rptr(i+1) - akeep%rptr(i))
      akeep%maxmn = max(akeep%maxmn, blkm, blkn)
      akeep%level(i) = akeep%level(akeep%sparent(i)) + 1
      info%maxfront = max(info%maxfront, blkn)
      info%maxdepth = max(info%maxdepth, akeep%level(i))
   end do

   ! Setup child_ptr, child_next and calculate work per subtree
   allocate(child_next(akeep%nnodes+1), child_head(akeep%nnodes+1), &
      akeep%child_ptr(akeep%nnodes+2), akeep%child_list(akeep%nnodes), &
      akeep%subtree_work(akeep%nnodes+1), stat=st)
   if (st .ne. 0) go to 490
   child_head(:) = -1
   do i = akeep%nnodes, 1, -1 ! backwards so child list is in order
      blkn = akeep%sptr(i+1) - akeep%sptr(i) 
      blkm = int(akeep%rptr(i+1) - akeep%rptr(i))
      j = akeep%sparent(i)
      ! Add to parent's child linked list
      child_next(i) = child_head(j)
      child_head(j) = i
      ! Calculate extra work at this node
      akeep%subtree_work(i) = 0
      do k = blkm, blkm-blkn+1, -1
         akeep%subtree_work(i) = akeep%subtree_work(i) + k**2
      end do
   end do
   akeep%subtree_work(akeep%nnodes+1) = 0
   ! Add work up tree, build child_ptr and child_list
   akeep%child_ptr(1) = 1
   do i = 1, akeep%nnodes+1
      if(i.lt.akeep%nnodes+1) then
         j = akeep%sparent(i)
         akeep%subtree_work(j) = akeep%subtree_work(j) + akeep%subtree_work(i)
      end if
      j = child_head(i)
      akeep%child_ptr(i+1) = akeep%child_ptr(i)
      do while(j.ne.-1)
         akeep%child_list(akeep%child_ptr(i+1)) = j
         akeep%child_ptr(i+1) = akeep%child_ptr(i+1) + 1
         j = child_next(j)
      end do
   end do

   ! Info
   info%matrix_rank = akeep%sptr(akeep%nnodes+1)-1
   info%num_sup = akeep%nnodes

   ! Store copy of info data in akeep
   akeep%flag = info%flag
   akeep%matrix_dup = info%matrix_dup
   akeep%matrix_missing_diag = info%matrix_missing_diag
   akeep%matrix_outrange = info%matrix_outrange
   akeep%maxdepth = info%maxdepth
   akeep%num_sup = info%num_sup
   akeep%ordering = info%ordering
   akeep%num_flops = info%num_flops

   return

   490 continue
   info%stat = st
   if (info%stat .ne. 0) then
      info%flag = MA97_ERROR_ALLOCATION
      call ma97_print_flag(context,nout,info%flag,st=info%stat)
   end if
 
end subroutine analyse_phase

!****************************************************************************
!
! Build a map from A to nodes
! lcol( nlist(2,i) ) = val( nlist(1,i) )
! nptr defines start of each node in nlist
!
subroutine build_map(n, ptr, row, perm, invp, nnodes, sptr, rptr, rlist, &
      nptr, nlist, st)
   ! Original matrix A
   integer, intent(in) :: n
   integer, dimension(n+1), intent(in) :: ptr
   integer, dimension(ptr(n+1)-1), intent(in) :: row
   ! Permutation and its inverse (some entries of perm may be negative to
   ! act as flags for 2x2 pivots, so need to use abs(perm))
   integer, dimension(n), intent(in) :: perm
   integer, dimension(n), intent(in) :: invp
   ! Supernode partition of L
   integer, intent(in) :: nnodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   ! Row indices of L
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   ! Output mapping
   integer, dimension(nnodes+1), intent(out) :: nptr
   integer, dimension(2, ptr(n+1)-1), intent(out) :: nlist
   ! Error check paramter
   integer, intent(out) :: st

   integer :: i, j, k, p
   integer(long) :: jj
   integer :: blkm
   integer :: col
   integer :: node
   integer, dimension(:), allocatable :: ptr2, row2, origin
   integer, dimension(:), allocatable :: map

   allocate(map(n), ptr2(n+3), row2(ptr(n+1)-1), origin(ptr(n+1)-1), stat=st)
   if(st.ne.0) return

   !
   ! Build transpose of A in ptr2, row2. Store original posn of entries in
   ! origin array.
   !
   ! Count number of entries in row i in ptr2(i+2). Don't include diagonals.
   ptr2(:) = 0
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         if (k.eq.i) cycle
         ptr2(k+2) = ptr2(k+2) + 1
      end do
   end do
   ! Work out row starts such that row i starts in posn ptr2(i+1)
   ptr2(1:2) = 1
   do i = 1, n
      ptr2(i+2) = ptr2(i+2) + ptr2(i+1)
   end do
   ! Drop entries into place
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         if (k.eq.i) cycle
         row2(ptr2(k+1)) = i
         origin(ptr2(k+1)) = j
         ptr2(k+1) = ptr2(k+1) + 1
      end do
   end do

   !
   ! Build nptr, nlist map
   !
   p = 1
   do node = 1, nnodes
      blkm = int(rptr(node+1) - rptr(node))
      nptr(node) = p

      ! Build map for node indices
      do jj = rptr(node), rptr(node+1)-1
         map(rlist(jj)) = int(jj-rptr(node)+1)
      end do

      ! Build nlist from A-lower transposed
      do j = sptr(node), sptr(node+1)-1
         col = invp(j)
         do i = ptr2(col), ptr2(col+1)-1
            k = abs(perm(row2(i))) ! row of L
            if (k<j) cycle
            nlist(2,p) = (j-sptr(node))*blkm + map(k)
            nlist(1,p) = origin(i)
            p = p + 1
         end do
      end do

      ! Build nlist from A-lower
      do j = sptr(node), sptr(node+1)-1
         col = invp(j)
         do i = ptr(col), ptr(col+1)-1
            k = abs(perm(row(i))) ! row of L
            if (k<j) cycle
            nlist(2,p) = (j-sptr(node))*blkm + map(k)
            nlist(1,p) = i
            p = p + 1
         end do
      end do
   end do
   nptr(nnodes+1) = p
   
end subroutine build_map

!****************************************************************************

!
! Factorize phase (no right-hand sides)
!
subroutine ma97_factor_double(matrix_type, val, akeep, fkeep, control, &
      info, scale, ptr, row)
   integer, intent(in) :: matrix_type 
      ! 3 for real, positive definite
      ! 4 for real, indefinite.
      !-3 for Hermitian, positive definite
      !-4 for Hermitian, indefinite.
      !-5 for complex, indefinite.
   real(wp), dimension(*), intent(in) :: val ! matrix values (lower triangle)
   type(ma97_akeep), intent(in)  :: akeep
   type(ma97_fkeep), intent(inout) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   real(wp), dimension(*), optional, intent(inout) :: scale ! used to hold
      ! row and column scaling factors. Must be set on entry if
      ! control%scaling <= 0
   integer, dimension(*), optional, intent(in) :: ptr ! must be present 
      ! if on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.
   integer, dimension(*), optional, intent(in) :: row ! must be present if
      ! on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.

   real(wp) :: x(1)

   call ma97_factor_solve_double(matrix_type, val, 0, x, 0, akeep, &
      fkeep, control, info, scale=scale, ptr=ptr, row=row)
end subroutine ma97_factor_double

!****************************************************************************
!
! Factorize phase plus solve for single x
!
subroutine ma97_factor_solve_one_double(matrix_type,val,x1,akeep, &
      fkeep,control,info,scale,ptr,row)
   integer, intent(in) :: matrix_type 
      ! 3 for real, positive definite
      ! 4 for real, indefinite.
      !-3 for Hermitian, positive definite
      !-4 for Hermitian, indefinite.
      !-5 for complex, indefinite.
   real(wp), dimension(*), intent(in) :: val ! matrix values (lower triangle)
   real(wp), intent(inout) :: x1(:) 
      ! must hold right-hand sides on entry. On exit, holds solution.
   ! For details of akeep, fkeep, control, info : see derived type description
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(inout) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   real(wp), dimension(*), optional, intent(inout) :: scale ! used to hold
      ! row and column scaling factors. Must be set on entry if
      ! control%scaling <= 0
   integer, dimension(*), optional, intent(in) :: ptr ! must be present
      ! if on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.
   integer, dimension(*), optional, intent(in) :: row ! must be present if
      ! on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.

   integer :: ldx1

   ldx1 = size(x1)

   call ma97_factor_solve_double(matrix_type,val,1,x1,ldx1,akeep, &
      fkeep,control,info,scale=scale,ptr=ptr,row=row)
end subroutine ma97_factor_solve_one_double

!****************************************************************************
!
! Factorize phase plus solve for multiple rhs
!
subroutine ma97_factor_solve_double(matrix_type, val, nrhs, x, ldx, &
      akeep, fkeep, control, info, scale, ptr, row)
   integer, intent(in) :: matrix_type 
      ! 3 for real, positive definite
      ! 4 for real, indefinite.
      !-3 for Hermitian, positive definite
      !-4 for Hermitian, indefinite.
      !-5 for complex, indefinite.
   real(wp), dimension(*), intent(in) :: val ! matrix values (lower triangle)
   integer :: nrhs ! number right-hand sides. must be at least 1
      ! if solve required
   integer :: ldx ! First extent of x. must be at least n.
   real(wp), intent(inout) :: x(ldx,nrhs) ! if nrhs > 0
      ! must hold right-hand sides on entry. On exit, holds solution.
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(inout) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   real(wp), dimension(akeep%n), optional, intent(inout) :: scale ! used to hold
      ! row and column scaling factors. Must be set on entry if
      ! control%scaling <= 0
   integer, dimension(akeep%n+1), optional, intent(in) :: ptr ! must be
      ! present if on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.
   integer, dimension(*), optional, intent(in) :: row ! must be present if
      ! on call to analyse phase, check = .false.. Must be unchanged 
      ! since that call.

   real(wp), dimension(:), allocatable :: val2
   character(len=50) :: context

   logical :: abort
   integer :: i
   logical :: lsolve ! set to true if want to solve at same time as factorize
   integer :: n, nz
   integer :: nout, nout1
   integer :: mp
   integer :: r
   logical :: sing
   integer :: st
   type(stack_type), dimension(:), target, allocatable :: stack ! one entry
      ! for each node, used to keep track of where memory is
   type(thread_stats), dimension(:), allocatable :: stats ! one copy
      ! per thread, accumulates per thread statistics that are then summed to
      ! obtain global stats in info.
   integer, dimension(:,:), allocatable :: map ! work array, one copy per
      ! thread. Size (0:n, num_threads), with 0 index used to track which
      ! node current map refers to.
   type(real_ptr_type), dimension(:), allocatable :: buf ! stores pointer
      ! to work array. One copy (and work array) per thread.
   integer :: num_threads, this_thread
   ! Solve parameters. Tree is broken up into multiple chunks. Parent-child
   ! relations between chunks are stored in fwd_ptr and fwd (see solve routine
   ! comments)
   integer :: nchunk, local_job
   integer, dimension(:), allocatable :: chunk_sa, chunk_en, fwd_ptr, fwd
   type(smalloc_type), pointer :: next_alloc

   lsolve = .true.
   if (ldx .eq. 0 .and. nrhs .eq. 0) lsolve = .false.

   ! Setup for any printing we may require
   context = 'ma97_factor'
   mp = control%unit_diagnostics
   if (lsolve) context = 'ma97_factor_solve'
   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   nout1 = control%unit_warning
   if (control%print_level < 0) nout1 = -1

   ! Perform appropriate printing
   if (control%print_level >= 1 .and. mp >= 0) then
      if (matrix_type == HSL_MATRIX_REAL_SYM_PSDEF) then
         write (mp,'(//a,i2,a)') &
            ' Entering ma97_factor with matrix_type = ', matrix_type, ' and :'
         write (mp,'(a,5(/a,i12),5(/a,es12.4))') &
            ' control parameters (control%) :', &
            ' print_level         Level of diagnostic printing           = ', &
            control%print_level,      &
            ' unit_diagnostics    Unit for diagnostics                   = ', &
            control%unit_diagnostics, &
            ' unit_error          Unit for errors                        = ', &
            control%unit_error,       &
            ' unit_warning        Unit for warnings                      = ', &
            control%unit_warning,     &
            ' scaling             Scaling control                        = ', &
            control%scaling
      else
         write (mp,'(//a,i2,a)') &
            ' Entering ma97_factor with matrix_type = ', matrix_type, ' and :'
         write (mp,'(a,i12,a,5(/a,i12),5(/a,es12.4))') &
            ' matrix_type                                                = ', &
              matrix_type,              &
            ' control parameters (control%) :', &
            ' print_level         Level of diagnostic printing           = ', &
            control%print_level,      &
            ' unit_diagnostics    Unit for diagnostics                   = ', &
            control%unit_diagnostics, &
            ' unit_error          Unit for errors                        = ', &
            control%unit_error,       &
            ' unit_warning        Unit for warnings                      = ', &
            control%unit_warning,     &
            ' scaling             Scaling control                        = ', &
            control%scaling,          &
            ' small               Small pivot size                       = ', &
            control%small,           &
            ' u                   Initial relative pivot tolerance       = ', &
            control%u,               &
            ' multiplier          Multiplier for increasing array sizes  = ', &
            control%multiplier
      end if

      if (lsolve .and. nrhs > 0) then
         write (mp,'(a,i12/a,i12)') &
            ' ldx                                                        = ', &
            ldx, &
            ' nrhs                                                       = ', &
            nrhs
      end if
   end if

   if (.not.allocated(akeep%sptr) .or. akeep%flag < 0) then
      ! Analyse cannot have been run
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      fkeep%flag = info%flag
      return
   end if

   ! Initialize info output
   info%flag = MA97_SUCCESS
   info%matrix_dup = akeep%matrix_dup
   info%matrix_missing_diag = akeep%matrix_missing_diag
   info%matrix_outrange = akeep%matrix_outrange
   info%maxdepth = akeep%maxdepth
   info%num_sup = akeep%num_sup
   info%ordering = akeep%ordering
   info%maxfront = 0
   info%num_neg = 0
   info%num_delay = 0
   info%num_factor = 0
   info%num_flops = 0
   info%num_sup = akeep%nnodes
   info%num_two = 0
   info%stat = 0

   fkeep%flag = 0
   st = 0

   n = akeep%n

   if (akeep%nnodes.eq.0) then
      info%flag = MA97_SUCCESS
      info%matrix_rank = 0
      fkeep%flag = info%flag
      return
   end if

   select case(matrix_type)
   case(HSL_MATRIX_REAL_SYM_PSDEF)
      fkeep%pos_def = .true.
   case(HSL_MATRIX_REAL_SYM_INDEF)
      fkeep%pos_def = .false.
   case default
      info%flag = MA97_ERROR_MATRIX_TYPE
      call MA97_print_flag(context,nout,info%flag)
      return
   end select

   ! If matrix has been checked, produce a clean version of val in val2
   if (akeep%check) then
      nz = akeep%ptr(n+1) - 1
      allocate (val2(nz),stat=st)
      if (st .ne. 0) go to 10
      call mc69_set_values(matrix_type, akeep%lmap, akeep%map, val, &
         nz, val2)
   else
      ! analyse run with no checking so must have ptr and row present
      if (.not.present(ptr)) info%flag = MA97_ERROR_PTR_ROW
      if (.not.present(row)) info%flag = MA97_ERROR_PTR_ROW
      if (info%flag < 0) then
         call ma97_print_flag(context,nout,info%flag)
         fkeep%flag = info%flag
         return
      end if
      nz = akeep%ne
   end if

   ! At this point, either  ptr, row, val   
   !                  or    akeep%ptr, akeep%row, val2
   ! hold the lower triangular part of A

   !
   ! Perform scaling if required
   !
   if (control%scaling.gt.0 .or. present(scale)) then
      if(allocated(fkeep%scaling)) then
         if(size(fkeep%scaling).lt.n) then
            deallocate(fkeep%scaling, stat=st)
            allocate(fkeep%scaling(n), stat=st)
         end if
      else
         allocate(fkeep%scaling(n), stat=st)
      end if
      if (st.ne.0) go to 10
   else
      deallocate(fkeep%scaling, stat=st)
   end if

   if(allocated(akeep%scaling) .and. control%scaling.lt.3) then
      info%flag = MA97_WARNING_MATCH_ORD_NO_SCALE
      call ma97_print_flag(context,nout1,info%flag)
   end if

   select case (control%scaling)
   case(:0) ! User supplied or none
      if (present(scale)) then
         do i = 1, n
            fkeep%scaling(i) = scale(akeep%invp(i))
         end do
      end if
   case(1) ! MC64
      if (akeep%check) then
         call mc64_scale(n, akeep%ptr, akeep%row, val2,  &
            fkeep%scaling, akeep%invp, st, control%action, sing)
      else
         call mc64_scale(n, ptr, row, val, fkeep%scaling, akeep%invp, &
            st, control%action, sing)
      end if
      if (st .ne. 0) go to 10

      if (sing) then
         info%flag = MA97_ERROR_SINGULAR
         call ma97_print_flag(context,nout,info%flag)
         fkeep%flag = info%flag
         return
      end if
      if (present(scale)) then
         do i = 1, n
            scale(akeep%invp(i)) = fkeep%scaling(i)
         end do
      end if

   case(2) ! MC77
      if (akeep%check) then
         call mc77_scale(n, akeep%invp, akeep%ptr, akeep%row, val2, &
            fkeep%scaling, info%flag77, st)
      else
         call mc77_scale(n, akeep%invp, ptr, row, val, fkeep%scaling, &
            info%flag77, st)
      end if
      if (st .ne. 0) go to 10
      if (info%flag77 < 0) then
         ! unexpected error
         info%flag = MA97_ERROR_MC77
         call ma97_print_flag(context,nout,info%flag)
         fkeep%flag = info%flag
         return
      end if
      if (present(scale)) then
         do i = 1, n
            scale(akeep%invp(i)) = fkeep%scaling(i)
         end do
      end if
   case(3) ! MC64 calculated during ordering
      if (.not.allocated(akeep%scaling)) then
         ! No scaling saved from analyse phase
         info%flag = MA97_ERROR_NO_SAVED_SCALING
         call ma97_print_flag(context,nout,info%flag)
         fkeep%flag = info%flag
         return
      end if
      do i = 1, n
         fkeep%scaling(i) = akeep%scaling(akeep%invp(i))
      end do
   case(4:) ! MC30
      if (akeep%check) then
         call mc30_scale(n, akeep%invp, akeep%ptr, akeep%row, val2, &
            fkeep%scaling, info%flag77, st)
      else
         call mc30_scale(n, akeep%invp, ptr, row, val, fkeep%scaling, &
            info%flag77, st)
      end if
      if (st .ne. 0) go to 10
      if (info%flag77 < 0) then
         ! unexpected error
         info%flag = MA97_ERROR_MC77
         call ma97_print_flag(context,nout,info%flag)
         fkeep%flag = info%flag
         return
      end if
      if (present(scale)) then
         do i = 1, n
            scale(akeep%invp(i)) = fkeep%scaling(i)
         end do
      end if
   end select

   ! Peform scaling of rhs if required
   if (allocated(fkeep%scaling) .and. nrhs.ge.1) then
      do r = 1, nrhs
         do i = 1, n
            x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
         end do
      end do
   end if
   !if(allocated(fkeep%scaling)) &
   !   print *, "minscale, maxscale = ", minval(fkeep%scaling), &
   !      maxval(fkeep%scaling)

   ! Setup data storage
   if(allocated(fkeep%nodes)) then
      if(size(fkeep%nodes).lt.akeep%nnodes+1) then
         deallocate(fkeep%nodes,stat=st)
         allocate(fkeep%nodes(akeep%nnodes+1), stat=st)
         if (st .ne. 0) go to 10
      end if
   else
      allocate(fkeep%nodes(akeep%nnodes+1), stat=st)
      if (st .ne. 0) go to 10
   end if
   fkeep%nodes(1:akeep%nnodes+1)%ndelay = 0
         
   allocate(stack(akeep%nnodes), stat=st)
   if (st .ne. 0) go to 10
   ! Do an inital storage allocation:
   ! * control%multiplier * n             integers (for nodes(:)%perm)
   ! * control%multiplier * (nfactor+2*n) reals    (for nodes(:)%lcol)
   if(associated(fkeep%alloc)) then
      if(fkeep%alloc%imem_size.lt. &
            max(n+0_long, int(control%multiplier*n,kind=long))) then
         deallocate(fkeep%alloc%imem, stat=st)
         fkeep%alloc%imem_size = &
            max(n+0_long, int(control%multiplier*n,kind=long))
         allocate(fkeep%alloc%imem(fkeep%alloc%imem_size),stat=st)
         if (st .ne. 0) go to 10
      end if
      if(fkeep%alloc%rmem_size.lt. max(akeep%nfactor+2*n, &
            int(control%multiplier*akeep%nfactor+2*n,kind=long))) then
         deallocate(fkeep%alloc%rmem, stat=st)
         fkeep%alloc%rmem_size = max(akeep%nfactor+2*n, &
            int(control%multiplier*akeep%nfactor+2*n,kind=long))
         allocate(fkeep%alloc%rmem(fkeep%alloc%rmem_size), stat=st)
         if (st .ne. 0) go to 10
      end if
      next_alloc => fkeep%alloc
      do while(associated(next_alloc))
         next_alloc%rhead = 0
         next_alloc%ihead = 0
         next_alloc => next_alloc%next_alloc
      end do
      nullify(fkeep%alloc%top_real, fkeep%alloc%top_int)

   else
      allocate(fkeep%alloc, stat=st)
      if (st .ne. 0) go to 10
      call smalloc_setup(fkeep%alloc, &
         max(n+0_long, int(control%multiplier*n,kind=long)), &
         max(akeep%nfactor+2*n, &
            int(control%multiplier*akeep%nfactor+2*n,kind=long)), st)
      if (st .ne. 0) go to 10
   end if

   num_threads = 1
!$ num_threads = omp_get_max_threads()

   ! If we only have a small amount of work then force serial execution
   if(akeep%num_flops.lt.control%factor_min) num_threads = 1

   allocate(stats(num_threads), map(0:n, num_threads), buf(num_threads), &
      stat=st)
   if (st .ne. 0) go to 10
   map(0, :) = -1 ! initally map unassociated with any node

   !%%if(dolog) open(unit=LOG_UNIT, file="ma97.log")
   !%%if(dograph) then
   !%%   open(unit=GRAPH_UNIT, file="ma97.dot")
   !%%   write(GRAPH_UNIT, "(a)") "graph tree {"
   !%%end if
   if(num_threads.gt.1) then
      ! Do actual factorization
      abort = .false.
!$OMP PARALLEL DEFAULT(NONE)                                                  &
!$OMP    PRIVATE(this_thread)                                                 &
!$OMP    SHARED(abort, akeep, buf, control, fkeep, i, info, ldx, map, n,   &
!$OMP       nrhs, stack, stats, val, val2, x)

      this_thread = 1
!$    this_thread = omp_get_thread_num()+1

      allocate(buf(this_thread)%val(akeep%maxmn*128), &
         stat=stats(this_thread)%st)
      if (stats(this_thread)%st .ne. 0) then
         abort = .true.
         stats(this_thread)%flag = MA97_ERROR_ALLOCATION
      end if
!$OMP SINGLE
      if (akeep%check) then
         if (allocated(fkeep%scaling)) then
            call inner_factor(1, akeep%nnodes+1, fkeep%pos_def, n,             &
               akeep%nnodes, akeep%child_ptr, akeep%child_list,                &
               akeep%subtree_work, akeep%invp, akeep%nptr, akeep%nlist,        &
               val2, fkeep%nodes, akeep%sptr, akeep%sparent, akeep%level,      &
               akeep%rptr, akeep%rlist, fkeep%alloc, stack, map, buf,          &
               nrhs, x, ldx, control, stats, abort, scale=fkeep%scaling)
         else
            call inner_factor(1, akeep%nnodes+1, fkeep%pos_def, n,             &
               akeep%nnodes, akeep%child_ptr, akeep%child_list,                &
               akeep%subtree_work, akeep%invp, akeep%nptr, akeep%nlist,        &
               val2, fkeep%nodes, akeep%sptr, akeep%sparent, akeep%level,      &
               akeep%rptr, akeep%rlist, fkeep%alloc, stack, map, buf,          &
               nrhs, x, ldx, control, stats, abort)
         end if
      else
         if (allocated(fkeep%scaling)) then
            call inner_factor(1, akeep%nnodes+1, fkeep%pos_def, n,             &
               akeep%nnodes, akeep%child_ptr, akeep%child_list,                &
               akeep%subtree_work, akeep%invp, akeep%nptr, akeep%nlist,        &
               val, fkeep%nodes, akeep%sptr, akeep%sparent, akeep%level,       &
               akeep%rptr, akeep%rlist, fkeep%alloc, stack, map, buf,          &
               nrhs, x, ldx, control, stats, abort, scale=fkeep%scaling)
         else
            call inner_factor(1, akeep%nnodes+1, fkeep%pos_def, n,             &
               akeep%nnodes, akeep%child_ptr, akeep%child_list,                &
               akeep%subtree_work, akeep%invp, akeep%nptr, akeep%nlist,        &
               val, fkeep%nodes, akeep%sptr, akeep%sparent, akeep%level,       &
               akeep%rptr, akeep%rlist, fkeep%alloc, stack, map, buf,          &
               nrhs, x, ldx, control, stats, abort)
         end if
      end if
      if (stats(this_thread)%st .ne. 0) then
         abort = .true.
         stats(this_thread)%flag = MA97_ERROR_ALLOCATION
      end if
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

   else ! running on a single thread

      allocate(buf(1)%val(akeep%maxmn*128), stat=st)
      if(st.ne.0) goto 10

      if (akeep%check) then
         if (allocated(fkeep%scaling)) then
            call subtree_factor(1, akeep%nnodes, fkeep%pos_def,               &
               akeep%child_ptr, akeep%child_list, n, akeep%invp, akeep%nptr,  &
               akeep%nlist, val2, akeep%nnodes, fkeep%nodes, akeep%sptr,      &
               akeep%sparent, akeep%level, akeep%rptr, akeep%rlist,           &
               fkeep%alloc, stack, map, 1, buf, nrhs, x, ldx, control,        &
               stats, 1, scale=fkeep%scaling)
         else
            call subtree_factor(1, akeep%nnodes, fkeep%pos_def,               &
               akeep%child_ptr, akeep%child_list, n, akeep%invp, akeep%nptr,  &
               akeep%nlist, val2, akeep%nnodes, fkeep%nodes, akeep%sptr,      &
               akeep%sparent, akeep%level, akeep%rptr, akeep%rlist,           &
               fkeep%alloc, stack, map, 1, buf, nrhs, x, ldx, control,        &
               stats, 1)
         end if
      else
         if (allocated(fkeep%scaling)) then
            call subtree_factor(1, akeep%nnodes, fkeep%pos_def,               &
               akeep%child_ptr, akeep%child_list, n, akeep%invp, akeep%nptr,  &
               akeep%nlist, val, akeep%nnodes, fkeep%nodes, akeep%sptr,       &
               akeep%sparent, akeep%level, akeep%rptr, akeep%rlist,           &
               fkeep%alloc, stack, map, 1, buf, nrhs, x, ldx, control,        &
               stats, 1, scale=fkeep%scaling)
         else
            call subtree_factor(1, akeep%nnodes, fkeep%pos_def,               &
               akeep%child_ptr, akeep%child_list, n, akeep%invp, akeep%nptr,  &
               akeep%nlist, val, akeep%nnodes, fkeep%nodes, akeep%sptr,       &
               akeep%sparent, akeep%level, akeep%rptr, akeep%rlist,           &
               fkeep%alloc, stack, map, 1, buf, nrhs, x, ldx, control,        &
               stats, 1)
         end if
      end if

   end if

   !%%if(dolog) close(LOG_UNIT)
   !%%if(dograph) then
   !%%   write(GRAPH_UNIT, "(a)") "}"
   !%%   close(GRAPH_UNIT)
   !%%end if

   ! Do reductions
   i = minval(stats(:)%flag)
   if(i.lt.0) then
      info%flag = i
      info%stat = maxval(stats(:)%st)
      if(info%stat.eq.0) info%stat = minval(stats(:)%st)
      st = info%stat
   end if
   i = max(info%flag, maxval(stats(:)%flag))
   info%maxfront = maxval(stats(:)%maxfront)
   info%num_factor = sum(stats(:)%num_factor)
   info%num_flops = sum(stats(:)%num_flops)
   info%num_delay = sum(stats(:)%num_delay)
   info%num_neg = sum(stats(:)%num_neg)
   info%num_two = sum(stats(:)%num_two)
   info%matrix_rank = akeep%sptr(akeep%nnodes+1)-1 - sum(stats(:)%num_zero)

   if (info%flag < 0) then
      call ma97_print_flag(context, nout, info%flag)
      fkeep%flag = info%flag
      return
   end if

   ! Do back subs if lsolve = .true.
   if (lsolve) then
      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
         write (control%unit_diagnostics,'(/a,i8,a)') &
            ' Factorisation complete. Solving for ',nrhs,' right-hand sides'
      ! We aim to have 4 chunks per thread to hopefully provide sufficient
      ! tree-level parallelism. Following call divides tree into these
      ! chunks.
      call slv_calc_chunk(akeep%nnodes, fkeep%nodes, akeep%sparent, &
         akeep%rptr, 4*num_threads, nchunk, chunk_sa, chunk_en, fwd_ptr, &
         fwd, st)
      if(st.ne.0) goto 10
      local_job = MA97_SOLVE_JOB_DIAG_BWD
      if(fkeep%pos_def) local_job = MA97_SOLVE_JOB_DIAG
      !call subtree_bwd_slv(akeep%nnodes, 1, local_job, fkeep%pos_def, &
      !   akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, &
      !   akeep%invp, nrhs, x, ldx, info%stat)
!$OMP PARALLEL DEFAULT(SHARED) IF(info%num_factor.ge.control%solve_min)
!$OMP SINGLE
      call bwd_slv_tasks(nchunk+1, chunk_sa, chunk_en, fwd_ptr, fwd, &
         local_job, fkeep%pos_def, akeep%nnodes, fkeep%nodes, &
         akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp, nrhs, x, ldx, &
         control%solve_blas3, st)
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
      if(st.ne.0) goto 10

      if (allocated(fkeep%scaling)) then
         do r = 1, nrhs
            do i = 1, n
               x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
            end do
         end do
      end if
   end if

   if(akeep%n.ne.info%matrix_rank) then
      ! Rank deficient
      ! Note: If we reach this point then must be control%action=.true.
      info%flag = MA97_WARNING_FACT_SINGULAR
      call ma97_print_flag(context, nout1, info%flag)
   end if

   if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      if (lsolve) then
         write (control%unit_diagnostics,'(/a)') &
            ' Completed factorisation and solve with:'
      else
         write (control%unit_diagnostics,'(/a)') &
            ' Completed factorisation with:'
      end if
      write (control%unit_diagnostics, &
         '(a,2(/a,i12),2(/a,es12.4),5(/a,i12))') &
         ' information parameters (info%) :', &
         ' flag                   Error flag                               = ',&
         info%flag, &
         ' maxfront               Maximum frontsize                        = ',&
         info%maxfront, &
         ' num_factor             Number of entries in L                   = ',&
         real(info%num_factor), &
         ' num_flops              Number of flops performed                = ',&
         real(info%num_flops), &
         ' num_two                Number of 2x2 pivots used                = ',&
         info%num_two, &
         ' num_delay              Number of delayed eliminations           = ',&
         info%num_delay, &
         ' rank                   Computed rank                            = ',&
         info%matrix_rank, &
         ' num_neg                Computed number of negative eigenvalues  = ',&
         info%num_neg
 
   end if

   fkeep%flag = info%flag
   fkeep%matrix_rank = info%matrix_rank
   fkeep%maxfront = info%maxfront
   fkeep%num_delay = info%num_delay
   fkeep%num_factor = info%num_factor
   fkeep%num_flops = info%num_flops
   fkeep%num_neg = info%num_neg
   fkeep%num_two = info%num_two
   return
   !!!!!!!!!!!!!!!!!!!!

   !
   ! Error handling
   !
   10 continue
   info%flag = MA97_ERROR_ALLOCATION
   info%stat = st
   fkeep%flag = info%flag
   call ma97_print_flag(context, nout, info%flag, st=info%stat)

end subroutine ma97_factor_solve_double

!**************************************************************
!
! Call mc64 to get a scaling, then symmetrize it
!
subroutine mc64_scale(n, ptr, row, val, scaling, invp, st, action, sing)
   integer, intent(in) :: n ! order of system
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   integer, dimension(n), intent(in) :: invp
   integer, intent(out) :: st ! stat parameter
   logical, intent(in) :: action ! controls action if matrix found to be
      ! singular
   logical, intent(out) :: sing ! set to true if matrix found to be singular

   type(mc64_control) :: control64
   type(mc64_info) :: info64

   integer :: i
   integer, dimension(:), allocatable :: perm64
   real(wp), dimension(:), allocatable :: cscale
   
   allocate(perm64(2*n), cscale(2*n), stat=st)
   if (st .ne. 0) return
   sing = .false.
   control64%checking = 1 ! checking disabled
   control64%ldiag = -1 ! printing disabled
   call mc64_matching(5, HSL_MATRIX_REAL_SYM_INDEF, n, n, ptr, row, val, &
      control64, info64, perm64, scale=cscale)
   select case(info64%flag)
   case(0)
      ! success; do nothing
   case(1)
      ! structually singular matrix
      if(.not.action) then
         ! abort
         sing = .true.
         return
      end if
   case(2)
      ! large scaling factors; do nothing
   case(-5)
      ! allocation error
      st = info64%stat
      return
   case default
      ! some other error; should never happen
      ! signal through stat=-99
      st = -99
      return
   end select

   do i = 1, n
      scaling(i) = exp(cscale(invp(i)))
   end do

end subroutine mc64_scale

!*******************************
!
! Following Ruiz and Ucar
! "A symmetry preserving algorithm for matrix scaling"
! we do one iteration of the inf norm and then 3 of the one norm.
!
subroutine mc77_scale(n, invp, ptr, row, val, scaling, info, st)
   integer, intent(in) :: n ! order of system
   integer, dimension(n), intent(in) :: invp
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: ne, i, j, k

   integer, allocatable, dimension(:) :: iw
   real(wp), allocatable, dimension(:) :: dw
   real(wp), dimension(:), allocatable :: val2

   integer :: icntl77(10), info77(10)
   real(wp) :: cntl77(10), rinfo77(10)

   ! Take absolute value of matrix to avoid overheads
   ne = ptr(n+1)-1
   allocate(val2(ne), iw(2*n), dw(2*n), stat=st)
   if (st .ne. 0) then
      info = MA97_ERROR_ALLOCATION
      return
   end if
   val2(1:ne) = abs(val(1:ne))

   ! Set controls
   call mc77id(icntl77, cntl77)
   icntl77(1) = -1 ! error messages
   icntl77(2) = -1 ! warning messages
   icntl77(3) = -1 ! diagnostic messages
   icntl77(4) = -1 ! disable checking
   icntl77(5) = 1 ! absolute value precomputed
   icntl77(6) = -1 ! symmetric matrix


   ! Single iteration of inf norm
   icntl77(7) = 1 ! max number of iterations
   call mc77ad(0, n, n, ne, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl77, cntl77, info77, rinfo77)

   info = info77(1)
   if (info < 0) return

   do i = 1, n
      scaling(i) = 1/dw(invp(i))
   end do

   ! Apply scaling
   ! Note: we could modify mc77 to take an input vector of scaling to make
   ! this step unnesscary
   do i = 1, n
      do j = ptr(i), ptr(i+1)-1
         k = row(j)
         val2(j) = scaling(i) * val2(j) * scaling(k)
      end do
   end do


   ! Up to 3 iterations of one norm
   icntl77(7) = 3 ! max number of iterations
   call mc77ad(1, n, n, ptr(n+1)-1, ptr, row, val2, iw, size(iw), &
      dw, size(dw), icntl77, cntl77, info77, rinfo77)

   info = info77(1)
   if (info < 0) return

   do i = 1, n
      scaling(i) = scaling(i) / dw(invp(i))
   end do
   
end subroutine mc77_scale


subroutine mc30_scale(n, invp, ptr, row, val, scaling, info, st)
   integer, intent(in) :: n ! order of system
   integer, dimension(n), intent(in) :: invp
   integer, intent(in) :: ptr(n+1) ! column pointers of A
   integer, intent(in) :: row(*) ! row indices of A (lower triangle)
   real(wp), intent(in) :: val(*) ! entries of A (in same order as in row).
   real(wp), dimension(n), intent(out) :: scaling
   integer, intent(out) :: info
   integer, intent(out) :: st

   integer :: nnz, i
   integer, allocatable, dimension(:) :: col
   real(wp), allocatable, dimension(:) :: work, scale30

   info = 0
   st = 0

   nnz = ptr(n+1)-1
   allocate(col(nnz), stat=st)
   if (st .ne. 0) then
      info = MA97_ERROR_ALLOCATION
      return
   end if

   do i = 1, n
      col(ptr(i):ptr(i+1)-1) = i
   end do

   allocate(work(4*n), scale30(n), stat=st)
   if (st .ne. 0) then
      info = MA97_ERROR_ALLOCATION
      return
   end if
   call mc30ad(n, nnz, val, row, col, scale30, work, 6, info)

   ! Exponentiate scaling factors
   do i = 1, n
      scaling(i) = exp(scale30(invp(i)))
   end do

   ! Avoid scaling if factors are too large
   if(maxval(scaling)>1d40) scaling(:) = 1
end subroutine mc30_scale

!****************************************************************************
!
! This subroutine performs the parallel decomposition at a tree level.
!
! For a given node it loops over the children and creates a task for each
! of them. Consecutive children are combined if they have less than
! control%min_subtree_work between them. If a child subtree contains more than
! control%min_subtree_work below it then a recursive call is made to this
! subroutine for the node at the top of that child subtree. Otherwise one or
! more children are passed to subtree_factor to do the actual main work.
! Once all children have been factored then subtree_factor is used to process
! this node only.
!
! This subroutine works recursively. The current subtree is defined by the
! nodes sa:en.
!
recursive subroutine inner_factor(sa, en, pos_def, n, nnodes, child_ptr, &
      child_list, subtree_work, invp, nptr, nlist, val, &
      nodes, sptr, sparent, level, rptr, rlist, alloc, stack, map, buf, &
      nrhs, x, ldx, control, stats, abort, scale, ntask, recursed)
   integer, intent(in) :: sa ! Start node of subtree
   integer, intent(in) :: en ! End node of subtree
   logical, intent(in) :: pos_def ! True if problem is supposedly pos-definite
   integer, intent(in) :: n
   integer, intent(in) :: nnodes
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer(long), dimension(*), intent(in) :: subtree_work ! Amount of work
      ! below a subtree
   integer, dimension(*), intent(in) :: invp
   integer, dimension(*), intent(in) :: nptr
   integer, dimension(2,*), intent(in) :: nlist
   real(wp), dimension(*), intent(in) :: val
   type(node_type), dimension(*), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: level
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   type(smalloc_type), target, intent(inout) :: alloc ! Contains actual
   ! memory allocated to contain L. Everything else points to it.
   type(stack_type), dimension(*), target :: stack ! Has an entry for every node
   integer, dimension(0:n,*), intent(out) :: map ! One copy for each thread
      ! work array containing mapping for node map(0) (which will change)
   type(real_ptr_type), dimension(*), intent(inout) :: buf ! Work
      ! array, one for each thread
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx, *), intent(inout) :: x
   type(ma97_control), intent(in) :: control
   type(thread_stats), dimension(*), intent(inout) :: stats ! thread-level
      ! info parameters, one for each thread.
   logical, intent(inout) :: abort ! Global abort flag.
   real(wp), dimension(*), optional, intent(in) :: scale
   integer, optional :: ntask ! Target number of tasks below this node
   integer, optional :: recursed ! Counts recursion depth to avoid stack
      ! overflows

   integer, parameter :: maxrecurse = 500 ! max levels of recursion to avoid
      ! stack overflow

   integer :: cnode
   integer :: i

   integer :: this_thread
   integer :: node

   !%%type(log_type) :: ltask

   integer :: current_start, current_end
   integer(long) :: current_work

   integer :: myntask, myrecursed
   integer :: num_threads

   !%%character(len=8) :: n1, n2

   if(abort) return

   num_threads = 1
!$ num_threads = omp_get_max_threads()

   ! Calculate number of tasks for subprocesses to use.
   ! If this is the root then set ntask to be 4 times the total available
   ! threads. Otherwise use whatever our parent told us to use.
   myntask = 4*num_threads
   if(present(ntask)) then
      myntask = ntask
   elseif(num_threads.eq.1) then
      myntask = 1
   end if

   if(present(recursed)) then
      myrecursed = recursed + 1
   else
      myrecursed = 1
   end if

   this_thread = 1
!$ this_thread = omp_get_thread_num()+1

   !print *, "INNER FACTOR", sa, en, "on", this_thread

   node = en

   !
   ! Loop over children. Combine them if together they have less than
   ! control%min_subtree_work, otherwise fire off last child as a task and
   ! start a new combined task with this child subtree.
   !

   current_start = sa
   current_end = sa-1
   current_work = 0
   do i = child_ptr(en), child_ptr(en+1)-1
      cnode = child_list(i)
      !print *, "   child", cnode
      ! Factorize subtree rooted at cnode
      if(current_work+subtree_work(cnode).gt.control%min_subtree_work) then
         ! Create a new task
         if(current_work.gt.0) then
            ! As we're going to have more than one task below this one,
            ! assign each to generate half as many tasks.
            if(current_start.eq.sa) myntask = myntask / 2
            !%%if(dograph) then
            !%%   call node_name(n1, en)
            !%%   call node_name(n2, current_end)
            !%%   write(GRAPH_UNIT, "(4a)") n1, " -- ", n2, ";"
            !%%end if
!$OMP TASK DEFAULT(NONE) IF(myntask.gt.1) &
!$OMP    FIRSTPRIVATE(current_end, current_start, current_work, myntask,      &
!$OMP       myrecursed)                                                       &
!%%!$OMP    PRIVATE(ltask)                                                    &
!$OMP    PRIVATE(this_thread)                                                 &
!$OMP    SHARED(abort, alloc, buf, child_list, child_ptr, control, invp,      &
!$OMP       ldx, level, map, n, nnodes, nodes, nlist, nptr, nrhs, num_threads,&
!$OMP       pos_def, rlist, rptr, scale, sparent, sptr, stack, stats,         &
!$OMP       subtree_work, x, val)
         this_thread = 1
!$       this_thread = omp_get_thread_num()+1
         if(.not.abort .and. (myrecursed.ge.maxrecurse .or. &
               current_work.le.control%min_subtree_work)) then
            !%%if(dolog .and. myntask.gt.1) &
            !%%   call log_start(ltask, "SF", current_end)
            call subtree_factor(current_start, current_end, pos_def,          &
               child_ptr, child_list, n, invp, nptr, nlist, val, nnodes,      &
               nodes, sptr, sparent, level, rptr, rlist, alloc, stack,        &
               map(:,this_thread), num_threads, buf, nrhs, x, ldx, control,   &
               stats, myntask, scale=scale)
            !%%if(dolog .and. myntask.gt.1) call log_stop(ltask, LOG_UNIT)
            if(stats(this_thread)%flag.lt.0) abort = .true.
            if(stats(this_thread)%st.ne.0) abort = .true.
         else
            call inner_factor(current_start, current_end, pos_def, n,         &
               nnodes, child_ptr, child_list, subtree_work, invp, nptr, nlist,&
               val, nodes, sptr, sparent, level, rptr, rlist, alloc, stack,   &
               map, buf, nrhs, x, ldx, control, stats, abort, scale=scale,    &
               ntask=myntask, recursed=myrecursed)
         end if
!$OMP END TASK
         end if
         current_start = current_end+1
         current_end = cnode
         current_work = subtree_work(cnode)
      else
         ! Merge with any previous outstanding children into a single task
         current_end = cnode
         current_work = current_work + subtree_work(cnode)
      end if
   end do
   if(current_end.ge.current_start) then
      ! Create a new task
      !%%if(dograph) then
      !%%   call node_name(n1, en)
      !%%   call node_name(n2, current_end)
      !%%   write(GRAPH_UNIT, "(4a)") n1, " -- ", n2, ";"
      !%%end if
!$OMP TASK DEFAULT(NONE) IF(myntask.gt.1) &
!$OMP    FIRSTPRIVATE(current_end, current_start, current_work, myntask,      &
!$OMP       myrecursed)                                                       &
!%%!$OMP    PRIVATE(ltask)                                                    &
!$OMP    PRIVATE(this_thread)                                                 &
!$OMP    SHARED(abort, alloc, buf, child_list, child_ptr, control, invp,      &
!$OMP       ldx, level, map, n, nnodes, nodes, nlist, nptr, nrhs, num_threads,&
!$OMP       pos_def, rlist, rptr, scale, sparent, sptr, stack, stats, val, x, &
!$OMP       subtree_work)
      this_thread = 1
!$    this_thread = omp_get_thread_num()+1
      if(.not. abort .and. ( myrecursed.ge.maxrecurse .or. &
            current_work.le.control%min_subtree_work)) then
         !%%if(dolog .and. myntask.gt.1) &
         !%%   call log_start(ltask, "SF", current_end)

         call subtree_factor(current_start, current_end, pos_def,             &
            child_ptr, child_list, n, invp, nptr, nlist, val, nnodes, nodes,  &
            sptr, sparent, level, rptr, rlist, alloc, stack,                  &
            map(:,this_thread), num_threads, buf, nrhs, x, ldx, control,      &
            stats, myntask, scale=scale)
         !%%if(dolog .and. myntask.gt.1) call log_stop(ltask, LOG_UNIT)
         if(stats(this_thread)%flag.lt.0) abort = .true.
         if(stats(this_thread)%st.ne.0) abort = .true.
      else
         call inner_factor(current_start, current_end, pos_def, n,            &
            nnodes, child_ptr, child_list, subtree_work, invp, nptr, nlist,   &
            val, nodes, sptr, sparent, level, rptr, rlist, alloc, stack, map, &
            buf, nrhs, x, ldx, control, stats, abort, scale=scale,            &
            ntask=myntask, recursed=myrecursed)
      end if
!$OMP END TASK
   end if

   ! Now wait for all child subtrees to complete before we try and assemble them
!$OMP TASKWAIT
   if(abort) return

   ! Reset number of threads for this level (if it has been changed)
   myntask = 2*num_threads
   if(present(ntask)) then
      myntask = ntask
   elseif(num_threads.eq.1) then
      myntask = 1
   end if

!   this_thread = 1
!!$ this_thread = omp_get_thread_num()+1

   if(node.le.nnodes) then
      !%%if(dolog) call log_start(ltask, "FS", node)
      call subtree_factor(node, node, pos_def, child_ptr, child_list, &
         n, invp, nptr, nlist, val, nnodes, nodes, sptr, sparent, level, &
         rptr, rlist, alloc, stack, map(:,this_thread), num_threads, buf, &
         nrhs, x, ldx, control, stats, &
         myntask, scale=scale)
      if(stats(this_thread)%flag.lt.0) abort = .true.
      if(stats(this_thread)%st.ne.0) abort = .true.
      !%%if(dolog) call log_stop(ltask, LOG_UNIT)
   end if
   return

end subroutine inner_factor

!*******************************
!
! This subroutine factorises the subtree(s) that include nodes sa through
! en. Any elements being passed to nodes numbered higher than en are allocated
! using Fortran allocate statemenets rather than stack_alloc.
!
! We maintain the factors seperately from the generated elements to avoid
! copying. Factors are stored in alloc, but pointed to by entries of nodes(:)
! for ease of reference.
!
! Generated elements are stored in a pair of stacks (need two so we can copy
! from child to parent). They are not touched until the factorization has
! been performed on columns that we expect to eliminate.
!
! Entries of A are only added to L just before they are expected to be
! eliminated.
!
subroutine subtree_factor(sa, en, pos_def, child_ptr, child_list, n, &
      invp, nptr, nlist, val, nnodes, nodes, sptr, sparent, level, rptr, &
      rlist, alloc, stack, map, num_threads, buf, nrhs, x, ldx, &
      control, stats, ntask, scale)
   integer, intent(in) :: sa ! Start node of subtree
   integer, intent(in) :: en ! End node of subtree
   logical, intent(in) :: pos_def ! True if problem is supposedly pos-definite
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   integer, dimension(*), intent(in) :: nptr
   integer, dimension(2,*), intent(in) :: nlist
   real(wp), dimension(*), intent(in) :: val
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(inout) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: level
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   type(smalloc_type), target, intent(inout) :: alloc ! Contains actual memory
      ! allocations for L. Everything else (within the subtree) is just a
      ! pointer to this.
   type(stack_type), dimension(*), target :: stack
   integer, dimension(0:*), intent(inout) :: map ! map refers to node specified
      ! by map(0), which is changed as we iterate. At task scheduling points
      ! this may change in non-obvious ways, so we need to track this carefully.
   ! explicit size required on buf to avoid gfortran-4.3 bug
   integer, intent(in) :: num_threads
   type(real_ptr_type), dimension(num_threads), intent(inout) :: buf
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx, *), intent(inout) :: x
   type(ma97_control), intent(in) :: control
   type(thread_stats), dimension(*), intent(inout) :: stats
   real(wp), dimension(*), optional, intent(in) :: scale
   integer, intent(in) :: ntask

   integer :: blkm
   integer :: blkn
   integer :: cblkm
   integer :: cblkn
   integer :: cn
   integer :: cm
   integer :: cnd
   integer :: cnode
   integer :: i
   integer(long) :: ii
   integer(long) :: ip
   integer(long) :: ip1
   integer :: j
   integer(long) :: jj
   integer :: k
   integer :: m
   integer :: nd
   integer :: ndelay
   integer :: nelim
   integer :: parent
   integer(long) :: stack_sz
   integer :: this_thread

   integer :: node
   real(wp), dimension(:), pointer :: lcol
   real(wp), dimension(:), pointer :: lstack
   integer, dimension(:), pointer :: lperm
   real(wp), dimension(:), allocatable :: xlocal
   type(stack_type), pointer :: lsptr
   type(stack_mem_type), pointer :: stack_mem_odd, stack_mem_even
   real(wp) :: dummy(1)

   this_thread = 1
!$ this_thread = omp_get_thread_num() + 1

   !print *, "SUBTREE ", sa, en, "on", omp_get_thread_num()

   nullify(stack_mem_odd, stack_mem_even)

   if(nrhs.ge.1) then
      allocate(xlocal(nrhs*n), stat=stats(this_thread)%st)
      if (stats(this_thread)%st.ne.0) go to 10
   end if

   do node = sa, en
      ! Setup basic data about node
      ndelay = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + ndelay
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      parent = sparent(node)
      !print *, "NODE ", node, "(", blkm, "x", blkn, ")"

      stats(this_thread)%maxfront = max(stats(this_thread)%maxfront, blkn)

      ! Allocate memory for L (lcol(1:blkm*blkn)), D (lcol(blkm*blkn+1:) and
      ! the local permutation (lperm).
      call smalloc(alloc, nodes(node)%lcol, (blkn+0_long)*blkm+2*blkn, &
         nodes(node)%rsmptr, nodes(node)%rsmsa, stats(this_thread)%st)
      if (stats(this_thread)%st .ne. 0) go to 10
      lcol => nodes(node)%lcol
      call smalloc(alloc, nodes(node)%perm, blkn+0_long, &
         nodes(node)%ismptr, nodes(node)%ismsa, stats(this_thread)%st)
      if (stats(this_thread)%st .ne. 0) go to 10
      lperm => nodes(node)%perm

      ! Zero node
      lcol(:) = 0

      ! Add A
      if (present(scale)) then
         do i = nptr(node), nptr(node+1)-1
            j = nlist(2,i)
            k = (j-1) / (blkm-ndelay) + 1 ! local column
            j = mod(j-1, blkm-ndelay) + 1 ! local row
            ip = (ndelay+k-1)*blkm + ndelay+j
            j = rlist( rptr(node)+j-1 )
            k = rlist( rptr(node)+k-1 )
            lcol(ip) = scale(j) * val(nlist(1,i)) * scale(k)
         end do
      else
         do i = nptr(node), nptr(node+1)-1
            j = nlist(2,i)
            k = (j-1) / (blkm-ndelay) ! column
            k = (ndelay+k)*blkm + mod(j-1, blkm-ndelay)+1 + ndelay
            lcol(k) = val(nlist(1,i))
         end do
      end if

      ! Build map only if we have children
      if (child_ptr(node).ne.child_ptr(node+1)) then
         ! Note: In parallel we may skip to a different node part way through
         ! using map, but not normally the case. Store which node map is
         ! associated with to be safe.
         map(0) = node
         do jj = rptr(node), rptr(node+1)-1
            map(rlist(jj)) = ndelay + int(jj-rptr(node))+1
         end do
      end if

      ! Add contribution from children
      nd = 0
      do cn = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(cn)
         lsptr => stack(cnode)
         !print *, "found child ", cnode
         cnd = nodes(cnode)%ndelay
         cblkn = sptr(cnode+1)-sptr(cnode)
         cblkm = int(rptr(cnode+1)-rptr(cnode)) + cnd
         m = blkm - blkn
         cm = cblkm - cblkn - cnd

         ! Multifrontal update - pull update off stack
         ! Loop over (expected) uneliminated columns of child
         if (associated(lsptr%val)) then
            ! Note: If all variables at child are delayed there is no update!
            do ii = rptr(cnode)+cblkn, rptr(cnode+1)-1
               if (map(rlist(ii)) .gt. blkn) cycle ! lstack gets done later
               ! Column goes into lcol; loop over rows
               ip = ( map(rlist(ii)) - 1 ) * blkm
               ip1 = ( ii - rptr(cnode)-cblkn ) * (cm+1) + 1
               do jj = ii, rptr(cnode+1)-1
                  k = map(rlist(jj))
                  lcol(ip+k) = lcol(ip+k) + lsptr%val(ip1)
                  ip1 = ip1 + 1
               end do
            end do
         end if

         ! Add delayed columns to map
         k = nd
         do i = nodes(cnode)%nelim+1, cblkn+cnd
            j = nodes(cnode)%perm(i)
            k = k + 1
            lperm(k) = j
            map( j ) = i-nodes(cnode)%nelim
         end do
         ip = nd*blkm
         ! Copy delayed columns into place
         do i = nodes(cnode)%nelim+1, cblkn+cnd
            ip1 = (i-1)*cblkm
            ! Delayed rows
            lcol(ip+nd+1:ip+nd+cblkn+cnd-nodes(cnode)%nelim) = &
               nodes(cnode)%lcol(ip1+nodes(cnode)%nelim+1:ip1+cblkn+cnd)
            ! Expected rows
            ip1 = ip1 + cblkn + cnd + 1
            do jj = rptr(cnode)+cblkn, rptr(cnode+1)-1
               k = map( rlist(jj) )
               lcol(ip+k) = nodes(cnode)%lcol(ip1)
               ip1 = ip1 + 1
            end do
            ip = ip + blkm
         end do
         nd = nd + cblkn+cnd - nodes(cnode)%nelim
      end do

      ! Initialise lperm
      j = ndelay + 1
      do i = sptr(node), sptr(node+1)-1
         lperm(j) = i
         j = j + 1
      end do

      !print *, "pre lcol = "
      !print "(5es12.4)", lcol
      !if (blkm.ge.10 .or. blkn.ge.10) &
      !   print *, "node sz ", blkm, blkn, "(", ndelay, ")"

      ! Factor the current node
      call rfact_block(pos_def, blkm, blkn, blkn, ndelay, &
         !lcol(1:blkm*blkn), &
         nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
         blkm, &
         !lcol(blkm*blkn+1:), &
         nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), &
         buf, &
         !lperm, &
         nodes(node)%ismptr%imem(nodes(node)%ismsa), &
         nelim, control, stats, &
         !lcol, &
         nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
         0, ntask)
      if(stats(this_thread)%flag.lt.0) return
      !print *, "Delayed", blkn-nelim

      do i = blkm, blkm-nelim+1, -1
         stats(this_thread)%num_factor = stats(this_thread)%num_factor + i
         stats(this_thread)%num_flops = stats(this_thread)%num_flops + i**2
      end do

      ! Record delays
      nodes(node)%nelim = nelim
      stats(this_thread)%num_delay = stats(this_thread)%num_delay + blkn - nelim

      !print *, "lperm = ", lperm
      !print *, "post = "
      !print "(4es12.4)", lcol

      ! Generate element if we are not a root (and there have been eliminations)
      if (blkm.ne.blkn) then
         stack_sz = (blkm-blkn+0_long)**2 + (blkm-blkn+0_long)*nrhs
         if(parent.le.en) then ! Is generated element in this subtree?
            lsptr => stack(node)
            ! Allocate space to store generated element
            ! To allow copying from child to parent but still operate using
            ! stacks, have different stacks for odd and even levels of the tree
            if (mod(level(node),2).eq.0) then
               call stack_alloc(stack_mem_even, lsptr%val, stack_sz, &
                  stats(this_thread)%st)
               if (stats(this_thread)%st .ne. 0) go to 10
               lsptr%stptr => stack_mem_even
               lsptr%stsa = stack_mem_even%head - stack_sz + 1
            else
               call stack_alloc(stack_mem_odd, lsptr%val, stack_sz, &
                  stats(this_thread)%st)
               if (stats(this_thread)%st .ne. 0) go to 10
               lsptr%stptr => stack_mem_odd
               lsptr%stsa = stack_mem_odd%head - stack_sz + 1
            end if
            lstack => lsptr%val
         else
            ! At root, need to pass generated element back to caller, so
            ! allocate using Fortran allocate
            lsptr => stack(node)
            allocate(lsptr%val(stack_sz), stat=stats(this_thread)%st)
            if (stats(this_thread)%st .ne. 0) go to 10
            lstack => lsptr%val
            nullify(lsptr%stptr)
         end if
      end if

      ! Calculate generated element
      ! G = buf^T lcol
      ! Note: This has the side effect of zeroing the element first
      m = blkm-blkn ! size of generated element
      if (m.gt.0 .and. nelim.gt.0) then
         ip = blkn + 1
         if(pos_def) then
            if(associated(lsptr%stptr)) then
               call mydsyrk(m, m, nelim, -one, &
                  !lcol(ip:), & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+ip-1), &
                  blkm, zero, &
                  !lstack, & ! equivalent to below line
                  lsptr%stptr%mem(lsptr%stsa), &
                  m, control%min_ldsrk_work)
            else
               call mydsyrk(m, m, nelim, -one, &
                  !lcol(ip:), & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+ip-1), &
                  blkm, zero, &
                  lstack, &
                  m, control%min_ldsrk_work)
            end if
         else
            if(associated(lsptr%stptr)) then
               call ldsrk(m, m, nelim, -one, &
                  !lcol(ip:), blkm, & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+ip-1), blkm, &
                  !lcol(blkm*blkn+1:), & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), &
                  zero, &
                  !lstack, m, &
                  lsptr%stptr%mem(lsptr%stsa), m, &
                  buf, control%min_ldsrk_work, stats)
            else
               call ldsrk(m, m, nelim, -one, &
                  !lcol(ip:), blkm, & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+ip-1), blkm, &
                  !lcol(blkm*blkn+1:), & ! equivalent to below line
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), &
                  zero, &
                  lstack, m, &
                  buf, control%min_ldsrk_work, stats)
            end if
            if(stats(this_thread)%flag.lt.0) go to 10
         end if
      elseif (m.gt.0) then
         lstack(:) = 0.0_wp
      end if

      if (child_ptr(node).ne.child_ptr(node+1) .and. map(0).ne.node) then
         ! Note: In parallel we may skip to a different node part way through
         ! using map. This is due to task scheduling points existing for node
         ! level parallelism.
         ! To get around this in the rare cases where it happens we store the
         ! node that the map is assoicated with in map(0) and restore it if
         ! it has gone wrong.
         ! Note: does not include the delays
         map(0) = node
         do jj = rptr(node), rptr(node+1)-1
            map(rlist(jj)) = ndelay + int(jj-rptr(node)+1)
         end do
      end if

      ! Perform forward solve(s) if required
      if(nrhs.gt.0) then
         ! Pull expected elimination variables from rhs vector
         do i = 0, nrhs-1
            do j = 1, blkn
               xlocal(j+i*blkn) = x(invp(lperm(j)), i+1)
            end do
         end do
         ! Update map to reflect new permutation from factorization
         do i = 1, blkn
            map(lperm(i)) = i
         end do
         ! Pull any updates to elimination variables from children
         do cn = child_ptr(node), child_ptr(node+1)-1
            cnode = child_list(cn)
            lsptr => stack(cnode)
            cnd = nodes(cnode)%ndelay
            cblkn = sptr(cnode+1)-sptr(cnode)
            cblkm = int(rptr(cnode+1)-rptr(cnode)) + cnd
            cm = cblkm - cblkn - cnd
            if (associated(lsptr%val)) then
               do ii = rptr(cnode)+cblkn, rptr(cnode+1)-1
                  if (map(rlist(ii)) .gt. blkn) cycle ! unelim'd vars done later
                  ip = map(rlist(ii))
                  ip1 = cm**2 + ii-(rptr(cnode)+cblkn)+1
                  do i = 0, nrhs-1
                     xlocal(ip+i*blkn) = &
                        xlocal(ip+i*blkn) + lsptr%val(ip1+i*cm)
                  end do
               end do
            end if
         end do

         ! Perform actual solve
         if(blkm-blkn.gt.0) then
            if(associated(stack(node)%stptr)) then
               call slv_fwd_mf(pos_def, nelim, &
                  !lcol, &
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
                  blkm, nrhs, xlocal, blkn, &
                  !lstack((blkm-blkn)**2+1:), & ! equivalent to below line
                  stack(node)%stptr%mem(stack(node)%stsa+(blkm-blkn)**2), &
                  blkm-blkn, control%solve_blas3)
            else
               call slv_fwd_mf(pos_def, nelim, &
                  !lcol, &
                  nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
                  blkm, nrhs, xlocal, blkn, &
                  lstack((blkm-blkn)**2+1:), & ! equivalent to below line
                  !stack(node)%val((blkm-blkn)**2+1:), &
                  blkm-blkn, control%solve_blas3)
            end if
         else
            call slv_fwd_mf(pos_def, nelim, lcol, blkm, nrhs, xlocal, blkn, &
               dummy, 0, control%solve_blas3)
         end if

         ! Copy eliminated (and delayed) variables back to x
         do j = 1, blkn
            x(invp(lperm(j)), 1:nrhs) = xlocal(j:j+(nrhs-1)*blkn:blkn)
         end do
      end if

      ! Add any contributions to generated element from children
      ! Note: We go backwards here to ensure freed in correct order
      do cn = child_ptr(node+1)-1, child_ptr(node), -1
         cnode = child_list(cn)
         lsptr => stack(cnode)
         !print *, "found child ", cnode
         cnd = nodes(cnode)%ndelay
         cblkn = sptr(cnode+1)-sptr(cnode)
         cblkm = int(rptr(cnode+1)-rptr(cnode)) + cnd
         ip1 = 1
         m = blkm - blkn
         cm = cblkm - cblkn - cnd

         ! Loop over (expected) uneliminated columns of child
         if (associated(lsptr%val)) then
            ! Note: If all variables at child are delayed there is no update!
            do ii = rptr(cnode)+cblkn, rptr(cnode+1)-1
               if (map(rlist(ii)) .le. blkn) cycle ! Already gone into lcol
               ! Column goes into lstack; loop over rows
               ip = ( map(rlist(ii)) - blkn - 1 ) * m
               ip1 = ( ii - rptr(cnode)-cblkn ) * (cm + 1) + 1
               do jj = ii, rptr(cnode+1)-1
                  k = map(rlist(jj)) - blkn
                  lstack(ip+k) = lstack(ip+k) + lsptr%val(ip1)
                  ip1 = ip1 + 1
               end do
               ! Add rhs contributions to stack too
               if(nrhs.gt.0) then
                  ip = m**2 + map(rlist(ii)) - blkn
                  ip1 = cm**2 + ii-(rptr(cnode)+cblkn)+1
                  do i = 0, nrhs-1
                     lstack(ip+i*m) = lstack(ip+i*m) + lsptr%val(ip1+i*cm)
                  end do
               end if
            end do
         end if

         ! Now remove old pointer
         if(associated(lsptr%stptr)) then
            stack_sz = size(lsptr%val) ! Explicit var to work around g95 bug
            if (mod(level(cnode),2).eq.0) then
               call stack_free(stack_mem_even, stack_sz);
            else
               call stack_free(stack_mem_odd, stack_sz);
            end if
         else
            deallocate(lsptr%val, stat=stats(this_thread)%st)
         end if
      end do

      ! Store generated element (and/or delays) in list for attention at parent
      if (blkm.ne.blkn) then
         if(parent.gt.en) then
            ! Need an atomic update as multiple tasks may try and update this
!$OMP ATOMIC
            nodes(parent)%ndelay = nodes(parent)%ndelay + blkn - nelim
         else
            ! Can be non-atomic as contained in subtree this task owns
            nodes(parent)%ndelay = nodes(parent)%ndelay + blkn - nelim
         end if
      end if

   end do
   !print *, "   done subtree ", sa, en, "on", omp_get_thread_num()

   return

   10 continue
   stats(this_thread)%flag = MA97_ERROR_ALLOCATION
   return
   
end subroutine subtree_factor

!*******************************
!
! This subroutine performs the forward solve corresponding to the subtree
! sa:en using a multifrontal style. Updates to variables eliminated at
! parent nodes are passed up the tree via a stack.
!
! This code is closely related to the subtree_factor routine above.
!
subroutine subtree_fwd_solve(sa, en, pos_def, child_ptr, child_list, n, &
      invp, nnodes, nodes, sptr, sparent, level, rptr, &
      rlist, stack, map, nrhs, x, ldx, blas3, st)
   integer, intent(in) :: sa
   integer, intent(in) :: en
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(nnodes), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: level
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   type(stack_type), dimension(*), target :: stack
   integer, dimension(0:*), intent(inout) :: map
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx, *), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st

   integer :: blkm
   integer :: blkn
   integer :: cblkm
   integer :: cblkn
   integer :: cn
   integer :: cm
   integer :: cnd
   integer :: cnode
   integer :: i
   integer(long) :: ii
   integer(long) :: ip
   integer(long) :: ip1
   integer :: j
   integer :: k
   integer(long) :: jj
   integer :: m
   integer :: ndelay
   integer :: nelim
   integer :: parent
   integer(long) :: stack_sz
   integer :: this_thread

   integer :: node
   real(wp), dimension(:), pointer :: lstack
   integer, dimension(:), pointer :: lperm
   real(wp), dimension(:), allocatable :: xlocal
   type(stack_type), pointer :: lsptr
   type(stack_mem_type), pointer :: stack_mem_odd, stack_mem_even
   real(wp) :: dummy(1)

   st = 0
   this_thread = 1
!$ this_thread = omp_get_thread_num() + 1

   !print *, "SUBTREE ", sa, en, "on", omp_get_thread_num()

   nullify(stack_mem_odd, stack_mem_even)

   if(nrhs.ge.1) then
      allocate(xlocal(nrhs*n), stat=st)
      if (st.ne.0) return
   end if

   do node = sa, en
      ! Setup basic data about node
      ndelay = nodes(node)%ndelay
      nelim = nodes(node)%nelim
      blkn = sptr(node+1) - sptr(node) + ndelay
      blkm = int(rptr(node+1) - rptr(node)) + ndelay
      parent = sparent(node)
      !print *, "NODE ", node, "(", blkm, "x", blkn, ")"

      lperm => nodes(node)%perm

      ! Build map only if we have children
      if (child_ptr(node).ne.child_ptr(node+1)) then
         ! Note: In parallel we may skip to a different node part way through
         ! using map, but not normally the case. Store which node map is
         ! associated with to be safe.
         map(0) = node
         do jj = rptr(node), rptr(node+1)-1
            map(rlist(jj)) = ndelay + int(jj-rptr(node))+1
         end do
      end if

      ! Generate element if we are not a root (and there have been eliminations)
      if (blkm.ne.blkn) then
         stack_sz = (blkm-blkn+0_long)*nrhs
         if(parent.le.en) then ! Is generated element in this subtree?
            lsptr => stack(node)
            ! Allocate space to store generated element
            ! To allow copying from child to parent but still operate using
            ! stacks, have different stacks for odd and even levels of the tree
            if (mod(level(node),2).eq.0) then
               call stack_alloc(stack_mem_even, lsptr%val, stack_sz, st)
               if (st .ne. 0) return
               lsptr%stptr => stack_mem_even
               lsptr%stsa = stack_mem_even%head - stack_sz + 1
            else
               call stack_alloc(stack_mem_odd, lsptr%val, stack_sz, st)
               if (st .ne. 0) return
               lsptr%stptr => stack_mem_odd
               lsptr%stsa = stack_mem_odd%head - stack_sz + 1
            end if
            lstack => lsptr%val
         else
            ! At root, need to pass generated element back to caller, so
            ! allocate using Fortran allocate
            lsptr => stack(node)
            allocate(lsptr%val(stack_sz), stat=st)
            if (st .ne. 0) return
            lstack => lsptr%val
            nullify(lsptr%stptr)
         end if
      end if

      !
      ! Perform forward solve(s) if required
      !
      ! Pull expected elimination variables from rhs vector
      do k = 1, nrhs
         ip = (k-1)*blkn+1
         do j = 1, blkn
            xlocal(ip) = x(invp(lperm(j)), k)
            ip = ip + 1
         end do
      end do
      ! Update map with new lperm
      do i = 1, blkn
         map(lperm(i)) = i
      end do
      ! Pull any updates to elimination variables from children
      do cn = child_ptr(node), child_ptr(node+1)-1
         cnode = child_list(cn)
         lsptr => stack(cnode)
         cnd = nodes(cnode)%ndelay
         cblkn = sptr(cnode+1)-sptr(cnode)
         cblkm = int(rptr(cnode+1)-rptr(cnode)) + cnd
         cm = cblkm - cblkn - cnd
         if (associated(lsptr%val)) then
            do k = 1, nrhs
               ip1 = (k-1)*cm
               do ii = rptr(cnode)+cblkn, rptr(cnode+1)-1
                  ip1 = ip1 + 1
                  if (map(rlist(ii)) .gt. blkn) cycle ! unelim'd vars done later
                  ip = map(rlist(ii)) + (k-1)*blkn
                  xlocal(ip) = xlocal(ip) + lsptr%val(ip1)
               end do
            end do
         end if
      end do

      ! Perform actual solve
      if(blkm-blkn.gt.0) then
         if(associated(stack(node)%stptr)) then
            call slv_fwd_mf(pos_def, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
               blkm, nrhs, xlocal, blkn, &
               !lstack, &
               stack(node)%stptr%mem(stack(node)%stsa), &
               blkm-blkn, blas3)
         else
            call slv_fwd_mf(pos_def, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
               blkm, nrhs, xlocal, blkn, &
               lstack, &
               !stack(node)%val, &
               blkm-blkn, blas3)
         end if
      else
         call slv_fwd_mf(pos_def, nelim, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), &
            blkm, nrhs, xlocal, blkn, &
            dummy, 0, blas3)
      end if

      ! Copy eliminated (and delayed) variables back to x
      do k = 1, nrhs
         ip = (k-1)*blkn
         do j = 1, blkn
            ip = ip + 1
            x(invp(lperm(j)), k) = xlocal(ip)
         end do
      end do

      ! Add any contributions to generated element from children
      ! Note: We go backwards here to ensure freed in correct order
      do cn = child_ptr(node+1)-1, child_ptr(node), -1
         cnode = child_list(cn)
         lsptr => stack(cnode)
         if (associated(lsptr%val)) then
            cnd = nodes(cnode)%ndelay
            cblkn = sptr(cnode+1)-sptr(cnode)
            cblkm = int(rptr(cnode+1)-rptr(cnode)) + cnd
            ip1 = 1
            m = blkm - blkn
            cm = cblkm - cblkn - cnd

            ! Loop over (expected) uneliminated columns of child
            ! Note: If all variables at child are delayed there is no update!
            do k = 1, nrhs
               ip1 = (k-1)*cm
               do ii = rptr(cnode)+cblkn, rptr(cnode+1)-1
                  ip1 = ip1 + 1
                  if (map(rlist(ii)) .le. blkn) cycle ! Already gone into lcol
                  ! Add rhs contributions to stack too
                  if(nrhs.gt.0) then
                     ip = map(rlist(ii)) - blkn + (k-1)*m
                     lstack(ip) = lstack(ip) + lsptr%val(ip1)
                  end if
               end do
            end do
         end if

         ! Now remove old pointer
         if(associated(lsptr%stptr)) then
            stack_sz = size(lsptr%val) ! explicit var to work around g95 bug
            if (mod(level(cnode),2).eq.0) then
               call stack_free(stack_mem_even, stack_sz);
            else
               call stack_free(stack_mem_odd, stack_sz);
            end if
         else
            deallocate(lsptr%val, stat=st)
         end if
      end do

   end do
   !print *, "   done subtree ", sa, en, "on", omp_get_thread_num()
end subroutine subtree_fwd_solve

!*************************************************************************
!
! This subroutine performs a forwards solve in a multifrontal fashion
! (i.e. a contribution is passed up the tree rather than written directly to
!  the right-hand-side vector)
!
subroutine slv_fwd_mf(pos_def, nelim, l, ldl, nrhs, x, ldx, xstack, ldxs, blas3)
   logical, intent(in) :: pos_def
   integer, intent(in) :: nelim
   real(wp), dimension(*), intent(in) :: l
   integer, intent(in) :: ldl
   integer, intent(in) :: nrhs
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: ldx
   real(wp), dimension(*), intent(out) :: xstack
   integer, intent(in) :: ldxs
   logical, intent(in) :: blas3

   ! Have a special case for nrhs=1 (Use BLAS2 rather than BLAS3)
   if(nrhs.eq.1 .and. .not.blas3) then
      ! Perform solve
      if(pos_def) then
         call dtrsv('Lower', 'Non-Trans', 'Non-Unit', nelim, l, ldl, x, 1)
      else
         call dtrsv('Lower', 'Non-Trans', 'Unit', nelim, l, ldl, x, 1)
      end if

      ! Update x values corresponding to delays
      if(ldx.gt.nelim) then
         call dgemv('N', ldx-nelim, nelim, -one, l(nelim+1), ldl, &
            x, 1, one, x(nelim+1), 1)
      end if

      ! Update uneliminated contributions
      if(ldxs.gt.0) then
         if(nelim.ne.0) then
            call dgemv('N', ldxs, nelim, -one, l(ldx+1), ldl, x, 1, &
               zero, xstack, 1)
         else
            xstack(1:ldxs) = zero
         end if
      end if
   else
      ! Perform solve
      if(pos_def) then
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Non-Unit', nelim, nrhs, &
            one, l, ldl, x, ldx)
      else
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Unit', nelim, nrhs, &
            one, l, ldl, x, ldx)
      end if

      ! Update x values corresponding to delays
      if(ldx.gt.nelim) then
         call dgemm('N', 'N', ldx-nelim, nrhs, nelim, -one, l(nelim+1), ldl, &
            x, ldx, one, x(nelim+1), ldx)
      end if

      ! Update uneliminated contributions
      if(ldxs.gt.0) then
         call dgemm('N', 'N', ldxs, nrhs, nelim, -one, l(ldx+1), ldl, x, ldx, &
            zero, xstack, ldxs)
      end if
   end if
end subroutine slv_fwd_mf

!*************************************************************************
!
! Following subroutine handles a symmetric L * L' multiplication by splitting
! it into tasks in order to expose parallelism.
! In the serial case could just use dsyrk (but don't for bit-compatibility).
! Split as:
! ( L_1 ) ( L_1' L_2' ) + ( C_11      )
! ( L_2 )                 ( C_21 C_22 )
! We avoid ops on C_12 as we are symmetric and don't care.
! Partition is by columns of C. [there may be more rows than columns]
! Recurse for ops on C_11 and C_22.
! Do ops on C_21 ourself.
!
! Note: n <= m
!
subroutine mydsyrk(m, n, k, alpha, l, ldl, beta, c, ldc, &
      min_ldsrk_work)
   integer, intent(in) :: m ! number of rows in c
   integer, intent(in) :: n ! number of cols in c
   integer, intent(in) :: k ! number of cols in l and ld
   real(wp), intent(in) :: alpha
   real(wp), dimension(*), intent(in) :: l
   integer, intent(in) :: ldl
   real(wp), intent(in) :: beta
   real(wp), dimension(*), intent(inout) :: c
   integer, intent(in) :: ldc
   integer, intent(in) :: min_ldsrk_work

   integer, parameter :: minsz = 64

   integer(long) :: ip
   !%%type(log_type) :: ltask
   integer :: blkn, i, blkm, j

   if(n.eq.0) return

   do i = 1, n, minsz
      blkn = min(minsz, n-i+1)
!$OMP TASK DEFAULT(NONE) IF(blkn*blkn*k.ge.min_ldsrk_work) &
!$OMP    FIRSTPRIVATE(blkn, i) &
!$OMP    PRIVATE(ip) &
!%%!$OMP       PRIVATE(ltask) &
!$OMP    SHARED(alpha, beta, c, k, l, ldc, ldl)
      !%%if(dolog) call log_start(ltask, "DL")
      ip = (i-1)*ldc + i
      call dsyrk('L', 'N', blkn, k, alpha, l(i), ldl, beta, c(ip), ldc)
      !%%if(dolog) call log_stop(ltask, LOG_UNIT)
!$OMP END TASK
      do j = i+blkn, m, minsz
         blkm = min(minsz, m-j+1)
!$OMP    TASK DEFAULT(NONE) IF(blkm*blkn*k.ge.min_ldsrk_work) &
!$OMP       FIRSTPRIVATE(i, j, blkn, blkm) &
!$OMP       PRIVATE(ip) &
!%%!$OMP       PRIVATE(ltask) &
!$OMP       SHARED(alpha, beta, c, k, l, ldc, ldl)
         !%%if(dolog) call log_start(ltask, "DL")
         ip = (i-1)*ldc + j
         call dgemm('N', 'T', blkm, blkn, k, alpha, l(j), ldl, &
           l(i), ldl, beta, c(ip), ldc)
         !%%if(dolog) call log_stop(ltask, LOG_UNIT)
!$OMP    END TASK
      end do
   end do
!$OMP TASKWAIT
end subroutine mydsyrk

!*************************************************************************
!
! Following subroutine handles a symmetric L * D * L' multiplication by
! splitting it into tasks to expose parallelism.
! Split as:
! ( L_1 ) ( L_1'D_1' L_2'D_2' ) + ( C_11      )
! ( L_2 )                         ( C_21 C_22 )
! We avoid ops on C_12 as we are symmetric and don't care.
! Partition is by columns of C. [there may be more rows than columns]
! Recurse for ops on C_11 and C_22.
! Do ops on C_21 ourself.
!
! Note: n <= m
!
subroutine ldsrk(m, n, k, alpha, l, ldl, d, beta, c, ldc, buf, &
      min_ldsrk_work, stats)
   integer, intent(in) :: m ! number of rows in c
   integer, intent(in) :: n ! number of cols in c
   integer, intent(in) :: k ! number of cols in l and ld
   real(wp), intent(in) :: alpha
   real(wp), dimension(*), target, intent(in) :: l
   integer, intent(in) :: ldl
   real(wp), dimension(*), intent(in) :: d
   real(wp), intent(in) :: beta
   real(wp), dimension(*), intent(inout) :: c
   integer, intent(in) :: ldc
   type(real_ptr_type), dimension(*), intent(inout) :: buf
   integer, intent(in) :: min_ldsrk_work
   type(thread_stats), dimension(*), intent(inout) :: stats

   integer, parameter :: minsz = 128!64

   integer(long) :: ip

   integer :: blkn, i, blkm, j
   integer :: this_thread
   !%%type(log_type) :: ltask

   if(n.eq.0) return

   do i = 1, n, minsz
      blkn = min(minsz, n-i+1)
      do j = i, m, minsz
         blkm = min(minsz, m-j+1)
!$OMP    TASK DEFAULT(NONE) IF(blkm*blkn*k.ge.min_ldsrk_work) &
!$OMP       FIRSTPRIVATE(i, j, blkn, blkm, alpha, beta, k, ldc, ldl, m, n) &
!$OMP       PRIVATE(ip, this_thread) &
!%%!$OMP       PRIVATE(ltask) &
!$OMP       SHARED(buf, c, d, l, stats)
         this_thread = 1
!$       this_thread = omp_get_thread_num()+1
         if(stats(this_thread)%flag.lt.0) goto 100
         !%%if(dolog) call log_start(ltask, "DL")
         ip = (i-1)*ldc + j
         ! Note that we may be switched to a thread that hasn't resized its
         ! buf array for this node yet.
         if(size(buf(this_thread)%val).lt.blkn*k) then
            deallocate(buf(this_thread)%val, stat=stats(this_thread)%st)
            allocate(buf(this_thread)%val(blkn*k), stat=stats(this_thread)%st)
            nullify(buf(this_thread)%chkptr)
            if(stats(this_thread)%st.ne.0) then
               stats(this_thread)%flag = MA97_ERROR_ALLOCATION
               goto 100
            end if
         end if
         if(.not.associated(buf(this_thread)%chkptr, target=l(i))) then
            call calc_ld(blkn, k, ldl, l(i), d, buf(this_thread)%val, &
               blkn)
            buf(this_thread)%chkptr => l(i)
         end if
         call dgemm('N', 'T', blkm, blkn, k, alpha, l(j), ldl, &
            buf(this_thread)%val, blkn, beta, c(ip), ldc)
         !%%if(dolog) call log_stop(ltask, LOG_UNIT)
         100 continue
!$OMP    END TASK
      end do
   end do
!$OMP TASKWAIT
end subroutine ldsrk

!*************************************************
!
! Subroutine takes a block of L and D^-1, calculates the
! corresponding block of LD
!
subroutine calc_ld(m, n, ldl, l_mat, d_vec, ld_mat, ldld)
   integer, intent(in) :: m ! number of rows in block
   integer, intent(in) :: n ! number of columns in block
   integer, intent(in) :: ldl ! leading dimension of block
   real(wp), dimension(n*ldl), intent(in) :: l_mat
   real(wp), dimension(2*n), intent(in) :: d_vec
   integer, intent(in) :: ldld
   real(wp), dimension(n*ldld), intent(out) :: ld_mat

   integer :: row, col
   integer :: lptr      ! index into l_mat
   integer :: ldptr     ! index into ld_mat
   real(wp) :: d1, d2, d3
   real(wp) :: det

   ! Loop over columns
   col = 1
   do while(col.le.n)
      lptr = 1 + (col-1)*ldl
      ldptr = 1 + (col-1)*ldld
      if(d_vec(2*col).eq.zero) then 
         ! 1x1 pivot
         if(d_vec(2*col-1).eq.zero) then
            d1 = zero
         else
            d1 = one/d_vec(2*col-1)
         end if
         do row = 1, m
            ld_mat(ldptr) = l_mat(lptr) * d1
            ldptr = ldptr + 1
            lptr  = lptr  + 1
         end do
         col = col   + 1
      else 
         ! 2x2 pivot
         d1 = d_vec(2*col - 1)
         d2 = d_vec(2*col)
         d3 = d_vec(2*col + 1)
         det = d1*d3 - d2*d2
         
         d1 = d1 / det
         d2 = d2 / det
         d3 = d3 / det
         ! Loop over rows
         do row = 1, m
            ld_mat(ldptr     ) =  d3*l_mat(lptr) - d2*l_mat(lptr+ldl)
            ld_mat(ldptr+ldld) = -d2*l_mat(lptr) + d1*l_mat(lptr+ldl)

            ldptr = ldptr + 1
            lptr  = lptr  + 1
         end do ! rows
         col   = col   + 2
      end if
   end do ! cols
end subroutine calc_ld

!*************************************************************************
!
! Solve phase single x. 
!
subroutine ma97_solve_one_double(x, akeep, fkeep, control, info, job)
   real(wp), dimension(:), intent(inout) :: x ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i) is the corresponding component of the
      ! right-hand side.On exit, if i has been used to index a variable,
      ! x(i) holds solution for variable i
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   integer, optional, intent(in) :: job

   integer :: ldx

   ldx = size(x)
   call ma97_solve_mult_double(1, x, ldx, akeep, fkeep, control, info, &
      job=job)
end subroutine ma97_solve_one_double

!*************************************************************************

subroutine ma97_solve_mult_double(nrhs, x, ldx, akeep, fkeep, control, &
      info, job)
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), intent(inout) :: x
   type(ma97_akeep), intent(in) :: akeep ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, if i has been used to index a variable,
      ! x(i,j) holds solution for variable i to system j
   ! For details of keep, control, info : see derived type description
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   integer, optional, intent(in) :: job ! used to indicate whether
      ! partial solution required
      ! job = 1 : forward eliminations only (PLX = B)
      ! job = 2 : diagonal solve (DX = B) (indefinite case only)
      ! job = 3 : backsubs only ((PL)^TX = B)
      ! job = 4 : diag and backsubs (D(PL)^TX = B) (indefinite case only)
      ! job absent: complete solve performed

   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i, r
   integer :: local_job ! local job parameter
   integer :: n
   integer :: nout
   type(stack_type), dimension(:), allocatable, target :: stack
   integer, dimension(:,:), allocatable :: map

   integer :: nchunk, num_threads
   integer, dimension(:), allocatable :: chunk_sa, chunk_en, fwd_ptr, fwd

   info%flag = MA97_SUCCESS

   ! Perform appropriate printing
   if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
         ' Entering ma97_solve with:'
      write (control%unit_diagnostics,'(a,4(/a,i12),(/a,i12))') &
         ' control parameters (control%) :', &
         ' print_level         Level of diagnostic printing        = ', &
         control%print_level, &
         ' unit_diagnostics    Unit for diagnostics                = ', &
         control%unit_diagnostics, &
         ' unit_error          Unit for errors                     = ', &
         control%unit_error, &
         ' unit_warning        Unit for warnings                   = ', &
         control%unit_warning, &
         ' nrhs                                                    = ', &
         nrhs
      if (nrhs > 1) write (control%unit_diagnostics,'(/a,i12)') &
         ' ldx                                                     = ', &
         ldx
   end if

   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   context = 'ma97_solve'

   if (akeep%nnodes.eq.0) return

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   info%flag = max(MA97_SUCCESS, fkeep%flag) ! Preserve warnings
   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   n = akeep%n
   if (ldx .lt. n) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' Increase ldx from ', ldx, ' to at least ', n
      return
   end if

   if (nrhs .lt. 1) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' nrhs must be at least 1. nrhs = ', nrhs
      return
   end if

   ! Copy previous phases' info data from akeep and fkeep
   info%matrix_dup = akeep%matrix_dup
   info%matrix_missing_diag = akeep%matrix_missing_diag
   info%matrix_outrange = akeep%matrix_outrange
   info%maxdepth = akeep%maxdepth
   info%matrix_rank = fkeep%matrix_rank
   info%maxfront = fkeep%maxfront
   info%num_delay = fkeep%num_delay
   info%num_factor = fkeep%num_factor
   info%num_flops = fkeep%num_flops
   info%num_neg = fkeep%num_neg
   info%num_sup = akeep%num_sup
   info%num_two = fkeep%num_two
   info%ordering = akeep%ordering

   ! Set local_job
   local_job = 0
   if (present(job)) then
      if (job .lt. 1 .or. job .gt. 4) info%flag = MA97_ERROR_JOB_OOR
      if (fkeep%pos_def .and. job == 2) info%flag = MA97_ERROR_JOB_OOR
      if (fkeep%pos_def .and. job == 4) info%flag = MA97_ERROR_JOB_OOR
      if (info%flag == MA97_ERROR_JOB_OOR) then
         call ma97_print_flag(context,nout,info%flag)
         return
      end if
      local_job = job
   end if

   if (allocated(fkeep%scaling)) then
      if (local_job == MA97_SOLVE_JOB_ALL .or. &
            local_job == MA97_SOLVE_JOB_FWD) then
         do r = 1, nrhs
            !x(1:n,r) = x(1:n,r) * fkeep%scaling(1:n)
            do i = 1, n
               x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
            end do
         end do
      end if
   end if

   ! We aim to have 4 chunks per thread to hopefully provide sufficient
   ! tree-level parallelism.
   num_threads = 1
!$ num_threads = omp_get_max_threads()
   call slv_calc_chunk(akeep%nnodes, fkeep%nodes, akeep%sparent, akeep%rptr, &
      4*num_threads, nchunk, chunk_sa, chunk_en, fwd_ptr, fwd, info%stat)
   if(info%stat.ne.0) goto 100

   if(control%solve_mf .and. ( local_job.eq.MA97_SOLVE_JOB_FWD .or. &
         local_job.eq.MA97_SOLVE_JOB_ALL) ) then
      ! use multifrontal forwards substitution (if doing supernodal then
      ! this gets done in inner_solve)
      allocate(map(0:akeep%n, num_threads), stack(akeep%nnodes), &
         stat=info%stat)
      if(info%stat.ne.0) goto 100
      !call subtree_fwd_solve(1, akeep%nnodes, fkeep%pos_def, &
      !   akeep%child_ptr, akeep%child_list, akeep%n, akeep%invp, &
      !   akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%sparent, &
      !   akeep%level, akeep%rptr, akeep%rlist, stack, map(:,1), nrhs, x, &
      !   ldx, info%stat)
!$OMP    PARALLEL DEFAULT(SHARED) IF(fkeep%num_factor.ge.control%solve_min)
!$OMP    SINGLE
      call fwd_slv_tasks(nchunk+1, chunk_sa, chunk_en, fwd_ptr, fwd, &
         fkeep%pos_def, akeep%child_ptr, akeep%child_list, akeep%n, &
         akeep%invp, akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%sparent, &
         akeep%level, akeep%rptr, akeep%rlist, stack, map, nrhs, x, ldx, &
         control%solve_blas3, info%stat)
!$OMP    END SINGLE NOWAIT
!$OMP    END PARALLEL
      if(info%stat.ne.0) goto 100
      ! Fudge local_job if required to perform backwards solve
      if(local_job.eq.MA97_SOLVE_JOB_ALL) then
         if(fkeep%pos_def) then
            local_job = MA97_SOLVE_JOB_BWD
         else
            local_job = MA97_SOLVE_JOB_DIAG_BWD
         end if
      elseif(local_job.eq.MA97_SOLVE_JOB_FWD) then
         local_job = -1 ! done
      end if
      !print *, "post fwd = ", x(1:3,1)
   end if

   ! Perform supernodal forward solve or diagonal solve (both in serial)
   call inner_solve(fkeep%pos_def, local_job, akeep%nnodes, &
      fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp, nrhs, &
      x, ldx, control%solve_blas3, info%stat)
   if (info%stat .ne. 0) goto 100

   if(local_job.eq.MA97_SOLVE_JOB_DIAG_BWD .or. &
         local_job.eq.MA97_SOLVE_JOB_BWD .or. &
         local_job.eq.MA97_SOLVE_JOB_ALL) then
      !call subtree_bwd_slv(akeep%nnodes, 1, local_job, fkeep%pos_def, &
      !   akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, &
      !   akeep%invp, nrhs, x, ldx, info%stat)
!$OMP PARALLEL DEFAULT(SHARED) IF(fkeep%num_factor.ge.control%solve_min)
!$OMP SINGLE
      call bwd_slv_tasks(nchunk+1, chunk_sa, chunk_en, fwd_ptr, fwd, &
         local_job, fkeep%pos_def, akeep%nnodes, fkeep%nodes, &
         akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp, nrhs, x, ldx, &
         control%solve_blas3, info%stat)
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
   end if
   if (info%stat .ne. 0) goto 100

   if (allocated(fkeep%scaling)) then
      if (local_job == MA97_SOLVE_JOB_ALL .or. &
            local_job == MA97_SOLVE_JOB_BWD .or. &
            local_job == MA97_SOLVE_JOB_DIAG_BWD) then
         do r = 1, nrhs
            !x(1:n,r) = x(1:n,r) * fkeep%scaling(1:n)
            do i = 1, n
               x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
            end do
         end do
      end if
   end if
   !print *, "soln = ", x(1:3,1)

   return

   100 continue
   info%flag = MA97_ERROR_ALLOCATION
   call ma97_print_flag(context,nout,info%flag,st=info%stat)
   return
end subroutine ma97_solve_mult_double

!*************************************************************************
!
! Alternative solve phase (A indefinite and singular).
!
subroutine ma97_solve_fredholm_double(nrhs, flag_out, x, ldx, &
      akeep, fkeep, control, info)
   integer, intent(in) :: nrhs
   logical, intent(out) :: flag_out(nrhs)
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,2*nrhs), intent(inout) :: x ! On entry, x must
      ! be set so that if i has been used to index a variable,
      ! x(i,j) is the corresponding component of the
      ! right-hand side for the jth system (j = 1,2,..., nrhs).
      ! On exit, x(1:n,1:nrhs) holds same solution as is computed by ma97_solve
      ! and, if flag_out*j)=.false., x(1:n,nrhs+j)
      ! holds the Fredholm alternative solution.
   ! For details of keep, control, info : see derived type description
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info

   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i, r
   integer :: n
   integer :: nout
   integer :: num_false ! number of inconsistent rhs
   integer :: num_rhs_bwd ! number of rhs for back sub (equal to nrhs if no
     ! no inconsistent rhs and is equal to 2*nrhs otherwise)
   type(stack_type), dimension(:), allocatable, target :: stack
   integer, dimension(:,:), allocatable :: map
   real(wp), dimension(:), allocatable :: tol
   real(wp), external :: dnrm2

   integer :: nchunk, num_threads
   integer, dimension(:), allocatable :: chunk_sa, chunk_en, fwd_ptr, fwd

   info%flag = MA97_SUCCESS
   ! initialise flag_out to true (consistent systems)
   flag_out(:) = .true.

   ! Perform appropriate printing
   if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
         ' Entering ma97_solve_fredholm with:'
      write (control%unit_diagnostics,'(a,4(/a,i12),(/a,i12))') &
         ' control parameters (control%) :', &
         ' print_level         Level of diagnostic printing        = ', &
         control%print_level, &
         ' unit_diagnostics    Unit for diagnostics                = ', &
         control%unit_diagnostics, &
         ' unit_error          Unit for errors                     = ', &
         control%unit_error, &
         ' unit_warning        Unit for warnings                   = ', &
         control%unit_warning, &
         ' nrhs                                                    = ', &
         nrhs
      if (nrhs > 1) write (control%unit_diagnostics,'(/a,i12)') &
         ' ldx                                                     = ', &
         ldx
   end if

   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   context = 'ma97_solve_fredholm'

   if (akeep%nnodes.eq.0) return

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   info%flag = max(MA97_SUCCESS, fkeep%flag) ! Preserve warnings
   ! immediate return if already had an error
   if (akeep%flag.lt.0 .or. fkeep%flag.lt.0) then
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (nrhs .lt. 1) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' nrhs must be at least 1. nrhs = ', nrhs
      return
   end if

   ! If matrix is non-singular, call ma97_solve
   if (fkeep%flag.ne.MA97_WARNING_FACT_SINGULAR) then
      call ma97_solve_mult_double(nrhs, x, ldx, akeep, fkeep, control, info)
      return
   end if

   n = akeep%n
   if (ldx .lt. n) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' Increase ldx from ', ldx, ' to at least ', n
      return
   end if

   ! Copy previous phases' info data from akeep and fkeep
   info%matrix_dup = akeep%matrix_dup
   info%matrix_missing_diag = akeep%matrix_missing_diag
   info%matrix_outrange = akeep%matrix_outrange
   info%maxdepth = akeep%maxdepth
   info%matrix_rank = fkeep%matrix_rank
   info%maxfront = fkeep%maxfront
   info%num_delay = fkeep%num_delay
   info%num_factor = fkeep%num_factor
   info%num_flops = fkeep%num_flops
   info%num_neg = fkeep%num_neg
   info%num_sup = akeep%num_sup
   info%num_two = fkeep%num_two
   info%ordering = akeep%ordering

   if (allocated(fkeep%scaling)) then
      do r = 1, nrhs
         do i = 1, n
            x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
         end do
      end do
   end if

    allocate(tol(nrhs),stat=info%stat)
    if(info%stat.ne.0) goto 100
    ! set tolerances for rhs's
    do i = 1,nrhs
       tol(i) = control%consist_tol * n * dnrm2(n,x(:,i),1)
    end do

   ! We aim to have 4 chunks per thread to hopefully provide sufficient
   ! tree-level parallelism.
   num_threads = 1
!$ num_threads = omp_get_max_threads()
   call slv_calc_chunk(akeep%nnodes, fkeep%nodes, akeep%sparent, akeep%rptr, &
      4*num_threads, nchunk, chunk_sa, chunk_en, fwd_ptr, fwd, info%stat)
   if(info%stat.ne.0) goto 100

   if(control%solve_mf) then
      ! use multifrontal forwards substitution (if doing supernodal then
      ! this gets done in inner_solve)
      allocate(map(0:akeep%n, num_threads), stack(akeep%nnodes), &
         stat=info%stat)
      if(info%stat.ne.0) goto 100

      !call subtree_fwd_solve(1, akeep%nnodes, fkeep%pos_def, &
      !   akeep%child_ptr, akeep%child_list, akeep%n, akeep%invp, &
      !   akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%sparent, &
      !   akeep%level, akeep%rptr, akeep%rlist, stack, map(:,1), nrhs, x, &
      !   ldx, info%stat)

!$OMP    PARALLEL DEFAULT(SHARED) IF(fkeep%num_factor.ge.control%solve_min)
!$OMP    SINGLE
      call fwd_slv_tasks(nchunk+1, chunk_sa, chunk_en, fwd_ptr, fwd, &
         fkeep%pos_def, akeep%child_ptr, akeep%child_list, akeep%n, &
         akeep%invp, akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%sparent, &
         akeep%level, akeep%rptr, akeep%rlist, stack, map, nrhs, x, ldx, &
         control%solve_blas3, info%stat)
!$OMP    END SINGLE NOWAIT
!$OMP    END PARALLEL
      if(info%stat.ne.0) goto 100

   else
      ! Perform supernodal forward solve (in serial)
      call inner_solve1(fkeep%pos_def, akeep%nnodes, &
         fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp, nrhs, &
         x, ldx, control%solve_blas3, info%stat)
      if (info%stat .ne. 0) goto 100
   end if

   ! Perform diagonal check for consistency and solve (in serial)
    
   call diag_solve(n, akeep%nnodes, fkeep%nodes, akeep%sptr, &
      akeep%rptr, akeep%invp, nrhs, x, ldx, tol, flag_out, num_false)

   ! Back solve is same in supernodal and multifrontal approach
   ! We do ONLY the back solve (no diagonal solve)

   !call subtree_bwd_slv(akeep%nnodes, 1, local_job, fkeep%pos_def, &
   !   akeep%nnodes, fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, &
   !   akeep%invp, nrhs, x, ldx, info%stat)

!$OMP PARALLEL DEFAULT(SHARED) IF(fkeep%num_factor.ge.control%solve_min)
!$OMP SINGLE
   if (num_false.eq.0) then
      ! no inconsistent right-hand sides so just do standard back sub.
      num_rhs_bwd = nrhs
   else
      ! back sub for both the consistent and inconsistent right-hand sides
      num_rhs_bwd = 2*nrhs
   end if
   call bwd_slv_tasks(nchunk+1, chunk_sa, chunk_en, fwd_ptr, fwd,    &
      MA97_SOLVE_JOB_BWD, fkeep%pos_def, akeep%nnodes,               &
      fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp,  &
      num_rhs_bwd, x, ldx, control%solve_blas3, info%stat)
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
   if (info%stat .ne. 0) goto 100

   if (allocated(fkeep%scaling)) then
      do r = 1, num_rhs_bwd
         do i = 1, n
            x(akeep%invp(i),r) = x(akeep%invp(i),r) * fkeep%scaling(i)
         end do
      end do
   end if

   return

   100 continue
   info%flag = MA97_ERROR_ALLOCATION
   call ma97_print_flag(context,nout,info%flag,st=info%stat)

   return

end subroutine ma97_solve_fredholm_double

!*************************************************************************
!
! Provide for multiplying a vector by L or L^T.
! y = Lx or y = L^Tx
!
subroutine ma97_lmultiply_one_double(trans, x1, y1, akeep, fkeep, &
      control, info)
   logical, intent(in) :: trans
   real(wp), dimension(:), intent(in) :: x1
   real(wp), dimension(:), intent(out) :: y1
   ! For details of keep, control, info : see derived type description
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info

   integer :: ldx, ldy ! needed for g95 with -i8 to avoid compiler bug

   ldx = size(x1)
   ldy = size(y1)
   call ma97_lmultiply_mult_double(trans,1,x1,ldx,y1,ldy, &
      akeep,fkeep,control,info)
end subroutine ma97_lmultiply_one_double

!*************************************************************************
!
! Provide for multiplying a vector by L or L^T.
! y = Lx or y = L^Tx
!
subroutine ma97_lmultiply_mult_double(trans, k, x, ldx, y, ldy, &
      akeep, fkeep, control, info)
   logical, intent(in) :: trans
   integer, intent(in) :: k
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,k), intent(in) :: x
   integer, intent(in) :: ldy
   real(wp), dimension(ldy,k), intent(out) :: y
   ! For details of keep, control, info : see derived type description
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info

   integer :: nout
   character(50)  :: context  ! Procedure name (used when printing).
   integer :: blkm, blkn, nd, nelim, node, i, r
   integer, dimension(:), allocatable :: map
   real(wp), dimension(:), allocatable :: xlocal
   real(wp), dimension(:,:), allocatable :: x2

   info%flag = MA97_SUCCESS
   info%stat = 0

   ! Perform appropriate printing
   if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
         ' Entering ma97_lmultiply with:'
      write (control%unit_diagnostics,'(a,5(/a,i12),(/a,l1))') &
         ' control parameters (control%) :', &
         ' print_level         Level of diagnostic printing        = ', &
         control%print_level, &
         ' unit_diagnostics    Unit for diagnostics                = ', &
         control%unit_diagnostics, &
         ' unit_error          Unit for errors                     = ', &
         control%unit_error, &
         ' unit_warning        Unit for warnings                   = ', &
         control%unit_warning, &
         ' k                                                       = ', &
         k, &
         ' trans                                                   = ', &
         trans
      if(k > 1) write (control%unit_diagnostics,'(2(/a,i12))') &
         ' ldx                                                     = ', &
         ldx, &
         ' ldy                                                     = ', &
         ldy
   end if

   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   context = 'ma97_lmultiply'

   if (akeep%nnodes.eq.0) return

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   info%flag = max(MA97_SUCCESS, fkeep%flag) ! Preserve warnings
   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (ldx.lt.akeep%n .or. ldy.lt.akeep%n) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (k .lt. 1) then
      info%flag = MA97_ERROR_X_SIZE
      call ma97_print_flag(context,nout,info%flag)
      if (nout .ge. 0) write (nout,'(a,i8,a,i8)') &
         ' k must be at least 1. k = ', k
      return
   end if

   allocate(xlocal((akeep%sptr(akeep%nnodes+1)-1)*k), &
      map(akeep%sptr(akeep%nnodes+1)-1), stat=info%stat)
   if(info%stat.ne.0) goto 100

   if(trans .and. allocated(fkeep%scaling)) then
      allocate(x2(akeep%n,k), stat=info%stat)
      if(info%stat.ne.0) goto 100
      do r = 1, k
         do i = 1, akeep%n
            x2(akeep%invp(i),r) = x(akeep%invp(i),r) / fkeep%scaling(i)
         end do
      end do
   end if

   y(1:akeep%n,:) = 0
   do node = 1, akeep%nnodes
      !print *, "node ", node
      nelim = fkeep%nodes(node)%nelim
      if (nelim.eq.0) cycle
      nd = fkeep%nodes(node)%ndelay
      blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
      blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
      !print *, "nelim, nd = ", nelim, nd
      !print *, "blkm, blkn = ", blkm, blkn
      
      if(allocated(x2)) then
         call node_mult(fkeep%pos_def, trans, k, x2, akeep%n, &
            y, ldy, akeep%rlist(akeep%rptr(node)), &
            akeep%invp, blkm, blkn, nelim, nd, &
            fkeep%nodes(node)%rsmptr%rmem(fkeep%nodes(node)%rsmsa), & ! lcol
            fkeep%nodes(node)%ismptr%imem(fkeep%nodes(node)%ismsa), & ! perm
            xlocal, map)
      else
         call node_mult(fkeep%pos_def, trans, k, x, ldx, y, ldy, &
            akeep%rlist(akeep%rptr(node)), &
            akeep%invp, blkm, blkn, nelim, nd, &
            fkeep%nodes(node)%rsmptr%rmem(fkeep%nodes(node)%rsmsa), & ! lcol
            fkeep%nodes(node)%ismptr%imem(fkeep%nodes(node)%ismsa), & ! perm
            xlocal, map)
      end if
   end do

   if(.not.trans .and. allocated(fkeep%scaling)) then
      do r = 1, k
         do i = 1, akeep%n
            y(akeep%invp(i),r) = y(akeep%invp(i),r) / fkeep%scaling(i)
         end do
      end do
   end if

   return

   100 continue
   info%flag = MA97_ERROR_ALLOCATION
   call ma97_print_flag(context,nout,info%flag,st=info%stat)
   return

end subroutine ma97_lmultiply_mult_double

!*************************************************************************

! Sparse forward solve (sparse rhs).
! Note: we pass in lflag and require x to be set to 0.0 on entry
! so we do not have loops of length n within the code (we need to 
! avoid loops of length n here since they will dominate total cost)

subroutine ma97_sparse_fwd_solve_double(nbi, bindex, b, order, lflag, &
      nxi, xindex, x, akeep, fkeep, control, info)
   integer, intent(in) :: nbi ! number of nonzero entries in right-hand side
   integer, intent(in) :: bindex(:) ! On entry, first nbi entries must 
      !  hold indices of  nonzero entries in the right-hand side. 
      !  Only first nbi entries are accessed. 
   real(wp), intent(in) :: b(:) ! If bindex(i)=k, b(k) must hold the k-th
      ! nonzero component of right-hand side; other entries of b not accessed.
   integer, intent(in) :: order(:) ! Must be passed unchanged form ma97_analyse.
   logical, intent(inout), dimension(:) :: lflag ! 
      ! Must be set by user on entry to false. used to flag which entries of x
      ! are non zero.
   integer, intent(out) :: nxi ! number of nonzero entries in the solution.
   integer, intent(out) :: xindex(:) ! First nxi entries holds indices of
      ! nonzero entries in solution. 
   real(wp), intent(inout) :: x(:) ! User must set to zero on entry.
      ! If xindex(i)=k, on exit x(k) holds the k-th
      ! nonzero component of solution; all other entries of x are zero.
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info

   character(50)  :: context  ! Procedure name (used when printing).
   integer :: i, j
   integer :: n
   integer :: nout

   info%flag = MA97_SUCCESS

   ! Perform appropriate printing
   if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
         ' Entering ma97_solve with:'
      write (control%unit_diagnostics,'(a,4(/a,i12),(/a,i12))') &
         ' control parameters (control%) :', &
         ' print_level         Level of diagnostic printing        = ', &
         control%print_level, &
         ' unit_diagnostics    Unit for diagnostics                = ', &
         control%unit_diagnostics, &
         ' unit_error          Unit for errors                     = ', &
         control%unit_error, &
         ' unit_warning        Unit for warnings                   = ', &
         control%unit_warning, &
         ' nbi                                                    = ', &
         nbi
   end if

   nout = control%unit_error
   if (control%print_level < 0) nout = -1
   context = 'ma97_sparse_fwd_solve'

   if (akeep%nnodes.eq.0) return

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   info%flag = max(MA97_SUCCESS, fkeep%flag) ! Preserve warnings
   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   n = akeep%n
   ! immediate return if nbi <= 0 .or. nbi > n
   if (nbi <= 0 .or. nbi > n) then
      info%flag = MA97_ERROR_NBI_OOR
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   ! Copy previous phases' info data from akeep and fkeep
   info%matrix_dup = akeep%matrix_dup
   info%matrix_missing_diag = akeep%matrix_missing_diag
   info%matrix_outrange = akeep%matrix_outrange
   info%maxdepth = akeep%maxdepth
   info%matrix_rank = fkeep%matrix_rank
   info%maxfront = fkeep%maxfront
   info%num_delay = fkeep%num_delay
   info%num_factor = fkeep%num_factor
   info%num_flops = fkeep%num_flops
   info%num_neg = fkeep%num_neg
   info%num_sup = akeep%num_sup
   info%num_two = fkeep%num_two
   info%ordering = akeep%ordering


   if (allocated(fkeep%scaling)) then
     do i = 1, nbi
        j = bindex(i)
        lflag(j) = .true.
        x(j) = b(j) * fkeep%scaling(order(j))
     end do
   else
      do i = 1, nbi
         j = bindex(i)
         lflag(j) = .true.
         x(j) = b(j)
      end do
   end if

   ! Perform supernodal forward solve (in serial)
   nxi = nbi
   xindex(1:nxi) = bindex(1:nbi)
   call inner_solve_sparse(fkeep%pos_def, akeep%nnodes, &
      fkeep%nodes, akeep%sptr, akeep%rptr, akeep%rlist, akeep%invp,  &
      x, n, nxi, xindex, control%solve_blas3, info%stat, lflag)

   if (info%stat .ne. 0) then
      info%flag = MA97_ERROR_ALLOCATION
      call ma97_print_flag(context,nout,info%flag,st=info%stat)
   end if

end subroutine ma97_sparse_fwd_solve_double

!*************************************************************************
!
! This subroutine provides a simple parallel wrapper around subtree_bwd_slv
! to spawn parallel tasks for different chunks of the assembly tree in a
! dependency-safe order.
!
recursive subroutine bwd_slv_tasks(chunk, chunk_sa, chunk_en, fwd_ptr, fwd, &
      job, pos_def, nnodes, nodes, sptr, rptr, rlist, invp, nrhs, x, &
      ldx, blas3, st)
   integer, intent(in) :: chunk
   integer, dimension(*), intent(in) :: chunk_sa
   integer, dimension(*), intent(in) :: chunk_en
   integer, dimension(*), intent(in) :: fwd_ptr
   integer, dimension(*), intent(in) :: fwd
   integer, intent(in) :: job ! controls whether we are doing forward
      ! eliminations, back substitutions etc.
   logical, intent(in) :: pos_def
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st  ! stat parameter

   integer :: i
   integer :: myst

   call subtree_bwd_slv(chunk_en(chunk), chunk_sa(chunk), job, pos_def, &
      nnodes, nodes, sptr, rptr, rlist, invp, nrhs, x, ldx, blas3, st)
   if(st.ne.0) return

   do i = fwd_ptr(chunk), fwd_ptr(chunk+1)-1
!$OMP TASK DEFAULT(NONE) &
!$OMP    FIRSTPRIVATE(i) &
!$OMP    PRIVATE(myst) &
!$OMP    SHARED(chunk_sa, chunk_en, fwd_ptr, fwd, job, pos_def, &
!$OMP       nnodes, nodes, sptr, rptr, rlist, invp, nrhs, x, ldx, blas3, st)
      if(st.eq.0) then
         call bwd_slv_tasks(fwd(i), chunk_sa, chunk_en, fwd_ptr, fwd, job, &
            pos_def, nnodes, nodes, sptr, rptr, rlist, invp, nrhs, x, &
            ldx, blas3, myst)
         if(myst.ne.0) st = myst
      end if
!$OMP END TASK
   end do

!$OMP TASKWAIT
   if(st.ne.0) return

end subroutine bwd_slv_tasks


!*************************************************************************
!
! This subroutine provides a simple parallel wrapper around subtree_fwd_slv
! to spawn parallel tasks for different chunks of the assembly tree in a
! dependency-safe order.
!
recursive subroutine fwd_slv_tasks(chunk, chunk_sa, chunk_en, fwd_ptr, fwd, &
      pos_def, child_ptr, child_list, n, invp, nnodes, nodes, sptr, sparent, &
      level, rptr, rlist, stack, map, nrhs, x, ldx, blas3, st)
   integer, intent(in) :: chunk
   integer, dimension(*), intent(in) :: chunk_sa
   integer, dimension(*), intent(in) :: chunk_en
   integer, dimension(*), intent(in) :: fwd_ptr
   integer, dimension(*), intent(in) :: fwd
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: child_ptr
   integer, dimension(*), intent(in) :: child_list
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   ! Note: gfortran-4.3 bug requires explicit size of nodes array
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sptr
   integer, dimension(*), intent(in) :: sparent
   integer, dimension(*), intent(in) :: level
   integer(long), dimension(*), intent(in) :: rptr
   integer, dimension(*), intent(in) :: rlist
   type(stack_type), dimension(*), target :: stack
   integer, dimension(0:n,*), intent(out) :: map
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx, *), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st

   integer :: i
   integer :: this_thread
   integer :: myst

   st = 0
   do i = fwd_ptr(chunk), fwd_ptr(chunk+1)-1
!$OMP TASK DEFAULT(NONE) &
!$OMP    FIRSTPRIVATE(i) &
!$OMP    PRIVATE(myst) &
!$OMP    SHARED(chunk_sa, chunk_en, fwd_ptr, fwd, pos_def, child_ptr, &
!$OMP       child_list, n, invp, nnodes, nodes, sptr, sparent, level, rptr, &
!$OMP       rlist, stack, map, nrhs, x, ldx, blas3, st)
      if(st.eq.0) then
         call fwd_slv_tasks(fwd(i), chunk_sa, chunk_en, fwd_ptr, fwd, &
            pos_def, child_ptr, child_list, n, invp, nnodes, nodes, sptr, &
            sparent, level, rptr, rlist, stack, map, nrhs, x, ldx, blas3, myst)
         if(myst.ne.0) st = myst
      end if
!$OMP END TASK
   end do

!$OMP TASKWAIT
   if(st.ne.0) return

   this_thread = 1
!$ this_thread = omp_get_thread_num() + 1

   call subtree_fwd_solve(chunk_sa(chunk), chunk_en(chunk), pos_def, &
      child_ptr, child_list, n, invp, nnodes, nodes, sptr, sparent, &
      level, rptr, rlist, stack, map(:,this_thread), nrhs, x, ldx, blas3, st)
   if(st.ne.0) return

end subroutine fwd_slv_tasks

!*************************************************************************
!
! This subroutine chops the assembly tree into chunks ready for parallel
! execution. The dependencies are encoded in fwd_ptr and fwd.
!
subroutine slv_calc_chunk(nnodes, nodes, sparent, rptr, ntask, nchunk, &
      chunk_sa, chunk_en, fwd_ptr, fwd, st)
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(*), intent(in) :: sparent
   integer(long), dimension(*), intent(in) :: rptr
   integer, intent(in) :: ntask ! Target number of chunks.
   integer, intent(out) :: nchunk ! Total number of chunks identified.
   integer, dimension(:), allocatable :: chunk_sa ! chunk_sa(i) is the first
      ! node in chunk i.
   integer, dimension(:), allocatable :: chunk_en ! chunk_en(i) is the last
      ! node in chunk i.
   integer, dimension(:), allocatable :: fwd_ptr ! children of chunk i are in
      ! fwd(fwd_ptr(i):fwd_ptr(i+1)-1)
   integer, dimension(:), allocatable :: fwd
   integer, intent(out) :: st

   integer, parameter :: min_entry = 10000 ! minimum number of entries in a
      ! chunk - if less than this merge into next consecutive chunk

   integer :: i, j
   integer :: nelim, nd, blkm
   integer :: parent
   integer :: node, chunk
   integer :: sa, en
   integer(long) :: num_factor, targ
   integer, dimension(:), allocatable :: first_child
   integer(long), dimension(:), allocatable :: ne
   integer, dimension(:), allocatable :: map
   integer, dimension(:), allocatable :: slv_parent

   ! Build first_child array and calculate global num_factor
   allocate(first_child(nnodes+1), ne(0:nnodes), map(nnodes+1), stat=st)
   if(st.ne.0) return
   first_child(:) = huge(first_child)
   num_factor = 0
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      nd = nodes(node)%ndelay
      blkm = int(rptr(node+1) - rptr(node)) + nd
      parent = sparent(node)

      num_factor = num_factor + nelim*blkm
      if(first_child(node).gt.node) first_child(node) = node
      first_child(parent) = min(first_child(node), first_child(parent))
   end do
   !print *, "Overall num_factor = ", num_factor

   ! Calculate target number of entries per chunk
   targ = num_factor / ntask
   !print *, "target = ", targ

   ! Iterate over nodes until current cumulative total exceeds target
   chunk = 1
   sa = 1
   num_factor = 0
   ne(0) = 0
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      nd = nodes(node)%ndelay
      blkm = int(rptr(node+1) - rptr(node)) + nd
      parent = sparent(node)

      num_factor = num_factor + nelim*blkm
      ne(node) = num_factor
      if(num_factor.ge.targ) then
         ! We have exceeded target number of entries. Split into disjoint
         ! subtrees using first children. Obey minimum size requirements,
         ! merging consecutive disjoint subtrees if necessary.
         i = node
         en = node
         do while(i.ge.sa)
            j = max(first_child(i), sa)
            if(ne(en) - ne(j-1) .ge. min_entry) then
               !print *, "Chunk ", j, en, "(dep ", sparent(en), ")", &
               !   ne(en) - ne(j-1)
               map(j:en) = chunk
               chunk = chunk + 1
               en = j-1
            end if
            i = j - 1
         end do
         if(en.ge.sa) then
            !print *, "Chunk ", sa, en, "(dep ", sparent(en), ")", &
            !   ne(en) - ne(j-1)
            map(sa:en) = chunk
            chunk = chunk + 1
         end if
         ne(node) = 0
         sa = node+1
         num_factor = 0
      end if
   end do
   if(sa.le.nnodes) then
      !print *, "Chunk ", sa, nnodes, " = ", num_factor
      map(sa:nnodes) = chunk
      chunk = chunk + 1
   end if
   map(nnodes+1) = chunk

   deallocate(ne, stat=st)
   deallocate(first_child, stat=st)
   
   ! We now need to manipulate data into an easy to use form in fwd_ptr, fwd
   nchunk = chunk-1
   allocate(chunk_sa(nchunk+1), chunk_en(chunk+1), fwd_ptr(nchunk+3), &
      fwd(nchunk+1), slv_parent(nchunk), stat=st)
   if(st.ne.0) return
   chunk_sa(nchunk+1) = -1; chunk_en(nchunk+1) = -2

   ! Set up chunk_sa, chunk_en and slv_parent (parent for chunk tree)
   ! Also count number of children using fwd_ptr at offset +2
   fwd_ptr(1:nchunk+3) = 0
   chunk = map(1)
   sa = 1
   do node = 1, nnodes+1
      if(map(node).eq.chunk) cycle
      ! Otherwise new chunk, dispatch old one
      !print *, "Chunk ", chunk, ":", sa, node-1, " depends on chunk ", &
      !   map(sparent(node-1))
      chunk_sa(chunk) = sa
      chunk_en(chunk) = node-1
      slv_parent(chunk) = map(sparent(node-1))
      fwd_ptr(slv_parent(chunk)+2) = fwd_ptr(slv_parent(chunk)+2) + 1
      sa = node
      chunk = map(node)
   end do

   fwd_ptr(1:2) = 1
   do chunk = 1, nchunk+1
      fwd_ptr(chunk+2) = fwd_ptr(chunk+2) + fwd_ptr(chunk+1)
   end do

   do chunk = 1, nchunk
      i = slv_parent(chunk)
      fwd(fwd_ptr(i+1)) = chunk
      fwd_ptr(i+1) = fwd_ptr(i+1) + 1
   end do

   !do chunk = 1, nchunk
   !   print *, "chunk ", chunk, ":", chunk_sa(chunk), chunk_en(chunk), &
   !      "slv_parent ", slv_parent(chunk)
   !   print *, "   fwd dep: ", fwd(fwd_ptr(chunk):fwd_ptr(chunk+1)-1)
   !end do
   !print *, "TOP LEVEL", fwd(fwd_ptr(nchunk+1):fwd_ptr(nchunk+2)-1)

end subroutine slv_calc_chunk

!*************************************************************************
!
! This subroutine performs a backwards solve on the chunk of nodes specified
! by sa:en.
!
subroutine subtree_bwd_slv(en, sa, job, pos_def, nnodes, nodes, sptr, &
      rptr, rlist, invp, nrhs, x, ldx, blas3, st)
   integer, intent(in) :: en
   integer, intent(in) :: sa
   logical, intent(in) :: pos_def
   integer, intent(in) :: job ! controls whether we are doing forward
      ! eliminations, back substitutions etc.
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st  ! stat parameter

   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: map

   st = 0
   allocate(xlocal(nrhs*(sptr(nnodes+1)-1)), map(sptr(nnodes+1)-1), stat=st)

!  print *, "SUBTREE ", en, sa, " on ", omp_get_thread_num()

   ! Backwards solve DL^Tx = z or L^Tx = z
   do node = en, sa, -1
      !print *, "node = ", node
      nelim = nodes(node)%nelim
      if (nelim.eq.0) cycle
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd
      
      if(nrhs.eq.1) then
         call slv_bwd_one(pos_def, job, rlist(rptr(node)), invp, x, &
            blkm, blkn, nelim, nd, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! node%lcol
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
            nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! node%perm
            xlocal, map, blas3)
      else
         call slv_bwd_mult(pos_def, job, rlist(rptr(node)), invp, &
            nrhs, x, ldx, blkm, blkn, nelim, nd, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! node%lcol
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
            nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! node%perm
            xlocal, map)
      end if
   end do

end subroutine subtree_bwd_slv

!*************************************************************************
!
! Provides serial versions of Forward (s/n) and diagonal solves.
!
subroutine inner_solve(pos_def, job, nnodes, nodes, sptr, rptr, rlist, &
      invp, nrhs, x, ldx, blas3, st)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job ! controls whether we are doing forward
      ! eliminations, back substitutions etc.
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st  ! stat parameter

   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: map

   st = 0
   allocate(xlocal(nrhs*(sptr(nnodes+1)-1)), map(sptr(nnodes+1)-1), stat=st)
   if(st.ne.0) return

   if (job == MA97_SOLVE_JOB_ALL .or. job == MA97_SOLVE_JOB_FWD) then
      ! Forwards solve Ly = b
      !print *, "pre = ", x
      do node = 1, nnodes
         nelim = nodes(node)%nelim
         if (nelim.eq.0) cycle
         nd = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + nd
         blkm = int(rptr(node+1) - rptr(node)) + nd
         
         if(nrhs.eq.1) then
            call slv_fwd_one(pos_def, rlist(rptr(node)), invp, x, &
               blkm, blkn, nelim, nd, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! nodes(node)%lcol
               nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! nodes(node)%perm
               xlocal, map, blas3)
         else
            call slv_fwd_mult(pos_def, rlist(rptr(node)), invp, nrhs, x, ldx, &
               blkm, blkn, nelim, nd, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! nodes(node)%lcol
               nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! nodes(node)%perm
               xlocal, map)
         end if
      end do
      !print *, "post fwd = ", x(1:3,1)
   end if

!   if (job.eq.MA97_SOLVE_JOB_ALL .or. &
!         job.eq.MA97_SOLVE_JOB_BWD .or. &
!         job.eq.MA97_SOLVE_JOB_DIAG_BWD) then
!
!      ! Backwards solve DL^Tx = z or L^Tx = z
!      do node = nnodes, 1, -1
!         !print *, "node = ", node
!         nelim = nodes(node)%nelim
!         if (nelim.eq.0) cycle
!         nd = nodes(node)%ndelay
!         blkn = sptr(node+1) - sptr(node) + nd
!         blkm = int(rptr(node+1) - rptr(node)) + nd
!         
!         if(nrhs.eq.1) then
!            call slv_bwd_one(pos_def, job, rlist(rptr(node)), invp, x, &
!               blkm, blkn, nelim, nd, &
!               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! node%lcol
!               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
!               nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! node%perm
!               xlocal, map, blas3)
!         else
!            call slv_bwd_mult(pos_def, job, rlist(rptr(node)), invp, &
!               nrhs, x, ldx, blkm, blkn, nelim, nd, &
!               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! node%lcol
!               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
!               nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! node%perm
!               xlocal, map)
!         end if
!      end do
!   end if

   if (job.eq.MA97_SOLVE_JOB_DIAG) then
      ! Diagonal solve Dx = z
      do node = nnodes, 1, -1
         nelim = nodes(node)%nelim
         if (nelim.eq.0) cycle
         nd = nodes(node)%ndelay
         blkn = sptr(node+1) - sptr(node) + nd
         blkm = int(rptr(node+1) - rptr(node)) + nd
         
         if(nrhs.eq.1) then
            call slv_diag_one(invp, x, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
               nodes(node)%ismptr%imem(nodes(node)%ismsa)) ! node%perm
         else
            call slv_diag_mult(invp, nrhs, x, ldx, nelim, &
               nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
               nodes(node)%ismptr%imem(nodes(node)%ismsa)) ! node%perm
         end if
      end do
      !print *, "post diag = ", x(1:3,1)
   end if

   !print *, "soln = ", x

end subroutine inner_solve

!*************************************************************************
!
! Provides serial version of Forward (supernodal).
!
subroutine inner_solve1(pos_def, nnodes, nodes, sptr, rptr, rlist, &
      invp, nrhs, x, ldx, blas3, st)
   logical, intent(in) :: pos_def
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,nrhs), intent(inout) :: x
   logical, intent(in) :: blas3
   integer, intent(out) :: st  ! stat parameter

   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: map

   st = 0
   allocate(xlocal(nrhs*(sptr(nnodes+1)-1)), map(sptr(nnodes+1)-1), stat=st)
   if(st.ne.0) return

   ! Forwards solve Ly = b
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      if (nelim.eq.0) cycle
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd
      
      if(nrhs.eq.1) then
         call slv_fwd_one(pos_def, rlist(rptr(node)), invp, x, &
            blkm, blkn, nelim, nd, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! nodes(node)%lcol
            nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! nodes(node)%perm
            xlocal, map, blas3)
      else
         call slv_fwd_mult(pos_def, rlist(rptr(node)), invp, nrhs, x, ldx, &
            blkm, blkn, nelim, nd, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! nodes(node)%lcol
            nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! nodes(node)%perm
            xlocal, map)
      end if
   end do
end subroutine inner_solve1

!*************************************************************************
!
! Sparse forward (s/n) solve.
!
subroutine inner_solve_sparse(pos_def, nnodes, nodes, sptr, rptr, rlist, &
      invp, x, n, nxi, xindex, blas3, st, lflag)
   logical, intent(in) :: pos_def
   integer, intent(in) :: nnodes
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(rptr(nnodes+1)-1), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: n
   real(wp), dimension(n), intent(inout) :: x
   integer, intent(inout) :: nxi
   integer, dimension(n), intent(out) :: xindex
   logical, intent(in) :: blas3
   integer, intent(out) :: st  ! stat parameter
   logical, dimension(n), intent(inout) :: lflag ! flags non zeros entries of x

   integer :: blkm
   integer :: blkn
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp), dimension(:), allocatable :: xlocal
   integer, dimension(:), allocatable :: map

   st = 0
   allocate(xlocal(sptr(nnodes+1)-1), map(sptr(nnodes+1)-1), stat=st)
   if(st.ne.0) return

   ! loop over nodos of tree.
   do node = 1, nnodes
      nelim = nodes(node)%nelim
      if (nelim.eq.0) cycle
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd        
      call slv_fwd_sparse(pos_def, rlist(rptr(node)), invp, x, &
            blkm, blkn, nelim, nd, &
            nodes(node)%rsmptr%rmem(nodes(node)%rsmsa), & ! nodes(node)%lcol
            nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! nodes(node)%perm
            xlocal, map, blas3, lflag, nxi, xindex)
   end do
end subroutine inner_solve_sparse

!*************************************************************************
!
! Provides serial version  of diagonal solve.
!
subroutine diag_solve(n, nnodes, nodes, sptr, rptr, invp, nrhs, x, ldx, &
      tol, flag_out, num_false)
   integer, intent(in) :: nnodes
   integer, intent(in) :: n
   type(node_type), dimension(*), intent(in) :: nodes
   integer, dimension(nnodes+1), intent(in) :: sptr
   integer(long), dimension(nnodes+1), intent(in) :: rptr
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,2*nrhs), intent(inout) :: x
   real(wp), intent (in) :: tol(nrhs) ! tolerance used to determine whether
       ! or not system is consistent
   logical, intent (inout) :: flag_out(nrhs)
   integer, intent (out) :: num_false ! number of inconsistent rhs

   integer :: blkm
   integer :: blkn
   integer :: i,j,r
   integer :: nd
   integer :: nelim
   integer :: node
   real(wp) :: temp

   num_false = 0
   do node = 1,nnodes
      nelim = nodes(node)%nelim
      if (nelim.eq.0) cycle
      nd = nodes(node)%ndelay
      blkn = sptr(node+1) - sptr(node) + nd
      blkm = int(rptr(node+1) - rptr(node)) + nd
      
      call check_diag(n, invp, nrhs, x, ldx, nelim, &
         nodes(node)%rsmptr%rmem(nodes(node)%rsmsa+blkm*blkn), & ! node%d
         nodes(node)%ismptr%imem(nodes(node)%ismsa), & ! node%perm
         flag_out,num_false,tol)

   end do

   if (num_false.lt.nrhs) then
     ! deal with any null rows and columns
     do i = sptr(nnodes+1),n
       j = invp(i) ! original row/col index
       do r = 1,nrhs
          if (.not. flag_out(r)) cycle
          temp = x(j,r)
          if (abs(temp) .gt. tol(r)) then
             flag_out(r) = .false.
             num_false = num_false + 1
             x(1:n,r+nrhs) = zero
             x(j,r+nrhs) = one
          end if
          if (num_false.eq.nrhs) exit
       end do
     end do
   end if

   if (num_false.gt.0 .and. num_false.lt.nrhs) then
     ! make sure that cols nrhs+1:2*nrhs are defined
     ! (as will perform back sub on these cols even if means
     ! we are solving for some consistent cols twice).
     do r = 1,nrhs
        if (.not. flag_out(r)) cycle
        x(1:n,r+nrhs) = 0
     end do
   end if

end subroutine diag_solve

!*************************************************************************
!
! Multiply single rhs by L
!
subroutine node_mult(pos_def, trans, nvec, x, ldx, y, ldy, rlist, invp, &
      blkm, blkn, nelim, nd, lcol, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   logical, intent(in) :: trans
   integer, intent(in) :: nvec
   real(wp), dimension(*), intent(in) :: x
   integer, intent(in) :: ldx
   real(wp), dimension(*), intent(inout) :: y
   integer, intent(in) :: ldy
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer :: i, k
   integer :: rp1
   character(len=1) :: diag
   
   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do
   
   ! Gather relevant part of x into xlocal
   if(.not.trans) then
      do k = 1, nvec
         do i = 1, nelim
            rp1 = map(i)
            xlocal(i+(k-1)*blkm) = x(rp1+(k-1)*ldx)
         end do
      end do
   else
      do k = 1, nvec
         do i = 1, blkm
            rp1 = map(i)
            xlocal(i+(k-1)*blkm) = x(rp1+(k-1)*ldx)
         end do
      end do
   end if

   diag = 'U'
   if(pos_def) diag = 'N'

   !print *, "lcol = ", lcol(1:blkm*blkn)
   !print *, "xlocal pre = ", xlocal(1:blkm)

   if(nvec.eq.1) then ! matrix-vector version
      if(.not.trans) then
         if (blkm-nelim.gt.0) &
            call dgemv('N', blkm-nelim, nelim, one, lcol(nelim+1), blkm, &
               xlocal, 1, zero, xlocal(nelim+1), 1)
         ! dtrmv comes last as it is destructive of xlocal
         call dtrmv('L', 'N', diag, nelim, lcol, blkm, xlocal, 1)
      else
         ! dtrmv comes first as dgemv is destructive of xlocal
         call dtrmv('L', 'T', diag, nelim, lcol, blkm, xlocal, 1)
         if (blkm-nelim.gt.0) &
            call dgemv('T', blkm-nelim, nelim, one, lcol(nelim+1), blkm, &
               xlocal(nelim+1), 1, one, xlocal, 1)
      end if
   else ! matrix-matrix version
      if(.not.trans) then
         if (blkm-nelim.gt.0) &
            call dgemm('N', 'N', blkm-nelim, nvec, nelim, one, lcol(nelim+1), &
               blkm, xlocal, blkm, zero, xlocal(nelim+1), blkm)
         ! dtrmv comes last as it is destructive of xlocal
         call dtrmm('Left', 'Lower', 'Non-T', diag, nelim, nvec, one, &
            lcol, blkm, xlocal, blkm)
      else
         ! dtrmv comes first as dgemv is destructive of xlocal
         call dtrmm('Left', 'Lower', 'Trans', diag, nelim, nvec, one, &
            lcol, blkm, xlocal, blkm)
         if (blkm-nelim.gt.0) &
            call dgemm('T', 'N', nelim, nvec, blkm-nelim, one, lcol(nelim+1), &
               blkm, xlocal(nelim+1), blkm, one, xlocal, blkm)
      end if
   end if

   !print *, "xlocal post = ", xlocal(1:blkm)

   ! Scatter solution back from xlocal
   if(.not.trans) then
      do k = 1, nvec
         do i = 1, blkm
            rp1 = map(i)
            y(rp1+(k-1)*ldy) = y(rp1+(k-1)*ldy) + xlocal(i+(k-1)*blkm)
         end do
      end do
   else
      do k = 1, nvec
         do i = 1, nelim
            rp1 = map(i)
            y(rp1+(k-1)*ldy) = xlocal(i+(k-1)*blkm)
         end do
      end do
   end if
end subroutine node_mult

!*************************************************************************
!
! Forward substitution single rhs
!
subroutine slv_fwd_one(pos_def, rlist, invp, x, blkm, blkn, nelim, nd, lcol, &
      lperm, xlocal, map, blas3)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map
   logical, intent(in) :: blas3

   integer(long) :: ip, ip2
   integer :: i, j, k
   integer :: rp1
   real(wp) :: ri, ri2
   
   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do
   
   ! Copy eliminated variables into xlocal
   do i = 1, nelim
      rp1 = map(i)
      xlocal(i) = x(rp1)
   end do

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Work with xlocal

      if(blas3) then
         if(pos_def) then
            call dtrsm('Left', 'Lower', 'Non-Trans', 'Non-Unit', nelim, 1, &
               one, lcol, blkm, xlocal, blkm)
         else
            call dtrsm('Left', 'Lower', 'Non-Trans', 'Unit', nelim, 1, &
               one, lcol, blkm, xlocal, blkm)
         end if
      else
         if(pos_def) then
            call dtrsv('L','N','N', nelim, lcol, blkm, xlocal, 1)
         else
            call dtrsv('L','N','U', nelim, lcol, blkm, xlocal, 1)
         end if
      end if

      if (blkm-nelim.gt.0) then
         if(blas3) then
            call dgemm('N', 'N', blkm-nelim, 1, nelim, -one, &
               lcol(nelim+1), blkm, xlocal, blkm, zero, &
               xlocal(nelim+1), blkm)
         else
            call dgemv('N', blkm-nelim, nelim, -one, lcol(nelim+1), blkm, &
               xlocal, 1, zero, xlocal(nelim+1), 1)
         end if
         ! Add contribution into x
         ! Delays first
         do i = nelim+1, blkm
            rp1 = map(i)
            x(rp1) = x(rp1) + xlocal(i)
         end do
      end if
   else
      do i = 1, nelim-1, 2
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         xlocal(i+1) = xlocal(i+1) - ri * lcol(ip+i+1)
         if(pos_def) xlocal(i+1) = xlocal(i+1) / lcol(ip2+i+1)
         ri2 = xlocal(i+1)
         do j = i+2, nelim
            xlocal(j) = xlocal(j) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
         do j = nelim+1, blkm
            rp1 = map(j)
            x(rp1) = x(rp1) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
      end do
      if(mod(nelim,2).eq.1) then
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         ! Commented loop redundant as i=nelim+1
         !do j = i+1, nelim
         !   xlocal(j) = xlocal(j) - ri * lcol(ip+j)
         !end do
         do j = nelim+1, blkm
            rp1 = map(j)
            x(rp1) = x(rp1) - ri * lcol(ip+j)
         end do
      end if
   end if

   ! Copy solution back from xlocal
   do i = 1, nelim
      rp1 = map(i)
      x(rp1) = xlocal(i)
   end do
end subroutine slv_fwd_one

!*************************************************************************
!
! Forward substitution multiple rhs
!
subroutine slv_fwd_mult(pos_def, rlist, invp, nrhs, x, ldx, blkm, blkn, nelim, &
      nd, lcol, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(blkm,*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip
   integer :: i, j, k, r
   integer :: rp1
   real(wp) :: ri
   
   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do
   
   ! Copy eliminated variables into xlocal
   do r = 1, nrhs
      do i = 1, nelim
         rp1 = map(i)
         xlocal(i,r) = x(rp1, r)
      end do
   end do

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Work with xlocal

      if(pos_def) then
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Non-Unit', nelim, nrhs, &
            one, lcol, blkm, xlocal, blkm)
      else
         call dtrsm('Left', 'Lower', 'Non-Trans', 'Unit', nelim, nrhs, &
            one, lcol, blkm, xlocal, blkm)
      end if

      if (blkm-nelim.gt.0) then
         call dgemm('N', 'N', blkm-nelim, nrhs, nelim, -one, &
            lcol(nelim+1), blkm, xlocal, blkm, zero, &
            xlocal(nelim+1,1), blkm)
         ! Add contribution into x
         do r = 1, nrhs
            ! Delays first
            do i = nelim+1, blkn
               rp1 = map(i)
               x(rp1,r) = x(rp1,r) + xlocal(i,r)
            end do
            ! Expected rows
            do j = blkn+1, blkm
               rp1 = map(j)
               x(rp1,r) = x(rp1,r) + xlocal(j,r)
            end do
         end do
      end if
   else
      do r = 1, nrhs
         do i = 1, nelim
            ip = (i-1)*blkm
            if(pos_def) xlocal(i,r) = xlocal(i,r) / lcol(ip+i)
            ri = xlocal(i,r)
            do j = i+1, nelim
               xlocal(j,r) = xlocal(j,r) - ri * lcol(ip+j)
            end do
            do j = nelim+1, blkm
               rp1 = map(j)
               x(rp1,r) = x(rp1,r) - ri * lcol(ip+j)
            end do
         end do
      end do
   end if

   ! Copy solution back from xlocal
   do r = 1, nrhs
      do i = 1, nelim
         rp1 = map(i)
         x(rp1,r) = xlocal(i,r)
      end do
   end do
end subroutine slv_fwd_mult

!*************************************************************************
!
! Forward substitution single SPARSE rhs
!
subroutine slv_fwd_sparse(pos_def, rlist, invp, x, blkm, blkn, nelim, nd, lcol,&
      lperm, xlocal, map, blas3, lflag, nxi, xindex)
   logical, intent(in) :: pos_def
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map
   logical, dimension(*), intent(inout) :: lflag ! lflag(j) = true if 
      !x(j) is non-zero
   logical, intent(in) :: blas3
   integer, intent(inout) :: nxi
   integer, dimension(*), intent(inout) :: xindex

   integer(long) :: ip, ip2
   integer :: i, j, k
   integer :: rp1
   real(wp) :: ri, ri2
   logical :: solve
   
   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1 + blkn - nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do
   
   ! Copy eliminated variables into xlocal. If corresponding entries of x
   ! are all zero, then nothing to do and return immediately.
   solve = .false.
   do i = 1,nelim
     rp1 = map(i)
     if (lflag(rp1)) solve = .true.
     xlocal(i) = x(rp1)
   end do
   if (.not. solve) return

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Work with xlocal

      if(blas3) then
         if(pos_def) then
            call dtrsm('Left', 'Lower', 'Non-Trans', 'Non-Unit', nelim, 1, &
               one, lcol, blkm, xlocal, blkm)
         else
            call dtrsm('Left', 'Lower', 'Non-Trans', 'Unit', nelim, 1, &
               one, lcol, blkm, xlocal, blkm)
         end if
      else
         if(pos_def) then
            call dtrsv('L','N','N', nelim, lcol, blkm, xlocal, 1)
         else
            call dtrsv('L','N','U', nelim, lcol, blkm, xlocal, 1)
         end if
      end if

      if (blkm-nelim.gt.0) then
         if(blas3) then
            call dgemm('N', 'N', blkm-nelim, 1, nelim, -one, &
               lcol(nelim+1), blkm, xlocal, blkm, zero, &
               xlocal(nelim+1), blkm)
         else
            call dgemv('N', blkm-nelim, nelim, -one, lcol(nelim+1), blkm, &
               xlocal, 1, zero, xlocal(nelim+1), 1)
         end if
         ! Add contribution into x
         ! Delays first
         do i = nelim+1, blkm
            rp1 = map(i)
            if (.not. lflag(rp1)) then
               lflag(rp1) = .true.
               nxi = nxi + 1
               xindex(nxi) = rp1
            end if
            x(rp1) = x(rp1) + xlocal(i)
         end do
      end if
   else
      do i = 1, nelim-1, 2
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         xlocal(i+1) = xlocal(i+1) - ri * lcol(ip+i+1)
         if(pos_def) xlocal(i+1) = xlocal(i+1) / lcol(ip2+i+1)
         ri2 = xlocal(i+1)
         do j = i+2, nelim
            xlocal(j) = xlocal(j) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
         do j = nelim+1, blkm
            rp1 = map(j)
            if (.not. lflag(rp1)) then
               lflag(rp1) = .true.
               nxi = nxi + 1
               xindex(nxi) = rp1
            end if
            x(rp1) = x(rp1) - ri * lcol(ip+j) - ri2 * lcol(ip2+j)
         end do
      end do
      if(mod(nelim,2).eq.1) then
         ip = (i-1)*blkm
         ip2 = i*blkm
         if(pos_def) xlocal(i) = xlocal(i) / lcol(ip+i)
         ri = xlocal(i)
         ! Commented loop redundant as i=nelim+1
         !do j = i+1, nelim
         !   xlocal(j) = xlocal(j) - ri * lcol(ip+j)
         !end do
         do j = nelim+1, blkm
            rp1 = map(j)
            if (.not. lflag(rp1)) then
               lflag(rp1) = .true.
               nxi = nxi + 1
               xindex(nxi) = rp1
            end if
            x(rp1) = x(rp1) - ri * lcol(ip+j)
         end do
      end if
   end if

   ! Copy solution back from xlocal
   do i = 1, nelim
      rp1 = map(i)
      if (.not. lflag(rp1)) then
         lflag(rp1) = .true.
         nxi = nxi + 1
         xindex(nxi) = rp1
      end if
      x(rp1) = xlocal(i)
   end do

end subroutine slv_fwd_sparse

!*************************************************************************
!
! Back substitution (with diagonal solve) single rhs
!
subroutine slv_bwd_one(pos_def, job, rlist, invp, x, blkm, blkn, nelim, &
      nd, lcol, d, lperm, xlocal, map, blas3)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job  ! used to indicate whether diag. sol. required
      ! job = 3 : backsubs only ((PL)^Tx = b)
      ! job = 0 or 4 : diag and backsubs (D(PL)^Tx = b)
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map
   logical, intent(in) :: blas3

   integer(long) :: ip
   integer :: i, j, k
   integer :: rp1, rp2

   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do

   if (job.eq.MA97_SOLVE_JOB_BWD .or. pos_def) then
      ! No diagonal solve. Copy eliminated variables into xlocal.
      do i = 1, nelim
         rp1 = map(i)
         xlocal(i) = x(rp1)
      end do
   else 
      ! Copy eliminated vars into xlocal while performing diagonal solve
      j = 1
      do while(j .le. nelim)
         if (d(2*j).ne.0) then
            ! 2x2 pivot
            rp1 = map(j)
            rp2 = map(j+1)
            xlocal(j)   = d(2*j-1) * x(rp1) + &
                          d(2*j)   * x(rp2)
            xlocal(j+1) = d(2*j)   * x(rp1) + &
                          d(2*j+1) * x(rp2)
            j = j + 2
         else
            ! 1x1 pivot
            if (d(2*j-1).eq.0.0_wp) then
               ! Zero pivot column
               xlocal(j) = 0.0_wp
            else
               ! Proper pivot
               rp1 = map(j)
               xlocal(j) = x(rp1) * d(2*j-1)
            end if
            j = j + 1
         end if
      end do
   end if
   !print *, "xlocal pre = ", xlocal(1:nelim) 

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      ! Delays
      do i = nelim+1, blkn
         rp1 = map(i)
         xlocal(i) = x(rp1)
      end do
      ! Expected rows
      do j = blkn+1, blkm
         rp1 = map(j)
         xlocal(j) = x(rp1)
      end do
      if (blkm-nelim.gt.0) then
         if(blas3) then
            call dgemm('Trans', 'Non-trans', nelim, 1, blkm-nelim, -one, &
               lcol(nelim+1), blkm, xlocal(nelim+1), blkm, one, xlocal, &
               blkm)
         else
            call dgemv('T', blkm-nelim, nelim, -one, lcol(nelim+1), blkm, &
               xlocal(nelim+1), 1, one, xlocal, 1)
         end if
      end if

      !print *, "xlocal2 pre = ", xlocal(1:nelim) 

      if(blas3) then
         if(pos_def) then
            call dtrsm('Left','Lower','Trans','Non-Unit', nelim, 1, one, &
               lcol, blkm, xlocal, blkm)
         else
            call dtrsm('Left','Lower','Trans','Unit', nelim, 1, one, lcol, &
               blkm, xlocal, blkm)
         end if
      else
         if(pos_def) then
            call dtrsv('L','T','N', nelim, lcol, blkm, xlocal, 1)
         else
            call dtrsv('L','T','U', nelim, lcol, blkm, xlocal, 1)
         end if
      end if

      !print *, "xlocal1 = ", xlocal(1:nelim)

      ! Copy solution back from xlocal
      do i = 1, nelim
         rp1 = map(i)
         x(rp1) = xlocal(i)
      end do
   else
      ! Do update with indirect addressing
      do i = 1, nelim
         ip = (i-1)*blkm
         do j = nelim+1, blkm
            rp1 = map(j)
            xlocal(i) = xlocal(i) - x(rp1) * lcol(ip+j)
         end do
      end do

      ! Solve with direct addressing
      if(pos_def) then
         do i = nelim, 1, -1
            ip = (i-1)*blkm
            rp1 = map(i)
            xlocal(i) = xlocal(i) - &
               sum(xlocal(i+1:nelim) * lcol(ip+i+1:ip+nelim))
            xlocal(i) = xlocal(i) / lcol(ip+i)
            x(rp1) = xlocal(i)
         end do
      else
         do i = nelim, 1, -1
            ip = (i-1)*blkm
            rp1 = map(i)
            xlocal(i) = xlocal(i) - &
               sum(xlocal(i+1:nelim) * lcol(ip+i+1:ip+nelim))
            x(rp1) = xlocal(i)
         end do
      end if
   end if
end subroutine slv_bwd_one

!*************************************************************************
!
! Back substitution (with diagonal solve) multiple rhs
!
subroutine slv_bwd_mult(pos_def, job, rlist, invp, nrhs, x, ldx, blkm, &
      blkn, nelim, nd, lcol, d, lperm, xlocal, map)
   logical, intent(in) :: pos_def
   integer, intent(in) :: job  ! used to indicate whether diag. sol. required
      ! job = 3 : backsubs only ((PL)^Tx = b)
      ! job = 0 or 4 : diag and backsubs (D(PL)^Tx = b)
   integer, dimension(*), intent(in) :: rlist
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: blkm
   integer, intent(in) :: blkn
   integer, intent(in) :: nelim
   integer, intent(in) :: nd
   real(wp), dimension(*), intent(in) :: lcol
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm
   real(wp), dimension(blkm,*), intent(out) :: xlocal
   integer, dimension(*), intent(out) :: map

   integer(long) :: ip
   integer :: i, j, k, r
   integer :: rp1, rp2

   do i = 1, blkn
      map(i) = invp( lperm(i) )
   end do
   k = 1+blkn-nd
   do i = blkn+1, blkm
      map(i) = invp( rlist(k) )
      k = k + 1
   end do

   if (job == MA97_SOLVE_JOB_BWD .or. pos_def) then 
      ! no diagonal solve. Copy eliminated variables into xlocal
      do r = 1, nrhs
         do i = 1, nelim
            rp1 = map(i)
            xlocal(i,r) = x(rp1, r)
         end do
      end do
   else 
      ! Copy eliminated vars into xlocal while performing diagonal solve
      do r = 1, nrhs
         j = 1
         do while(j .le. nelim)
            if (d(2*j).ne.0) then
               ! 2x2 pivot
               rp1 = map(j)
               rp2 = map(j+1)
               xlocal(j,r)   = d(2*j-1) * x(rp1,r) + &
                               d(2*j)   * x(rp2,r)
               xlocal(j+1,r) = d(2*j)   * x(rp1,r) + &
                               d(2*j+1) * x(rp2,r)
               j = j + 2
            else
               ! 1x1 pivot
               if (d(2*j-1).eq.0.0_wp) then
                  ! Zero pivot column
                  xlocal(j,r) = 0.0_wp
               else
                  ! Proper pivot
                  rp1 = map(j)
                  xlocal(j,r) = x(rp1,r) * d(2*j-1)
               end if
               j = j + 1
            end if
         end do
      end do
   end if

   ! Perform the solve
   if (blkm.gt.10 .and. nelim.gt.4) then
      do r = 1, nrhs
         ! Delays
         do i = nelim+1, blkn
            rp1 = map(i)
            xlocal(i,r) = x(rp1,r)
         end do
         ! Expected rows
         do j = blkn+1, blkm
            rp1 = map(j)
            xlocal(j,r) = x(rp1,r)
         end do
      end do
      if (blkm-nelim.gt.0) then
         call dgemm('Trans', 'Non-trans', nelim, nrhs, blkm-nelim, -one, &
            lcol(nelim+1), blkm, xlocal(nelim+1,1), blkm, one, xlocal, &
            blkm)
      end if

      if(pos_def) then
         call dtrsm('Left','Lower','Trans','Non-Unit', nelim, nrhs, one, lcol, &
            blkm, xlocal, blkm)
      else
         call dtrsm('Left','Lower','Trans','Unit', nelim, nrhs, one, lcol, &
            blkm, xlocal, blkm)
      end if
      do r = 1, nrhs
         do i = 1, nelim
            rp1 = map(i)
            x(rp1,r) = xlocal(i,r)
         end do
      end do
   else
      ! Do update with indirect addressing
      do r = 1, nrhs
         do i = 1, nelim
            ip = (i-1)*blkm
            do j = nelim+1, blkm
               rp1 = map(j)
               xlocal(i,r) = xlocal(i,r) - x(rp1,r) * lcol(ip+j)
            end do
         end do

         ! Solve with direct addressing
         if(pos_def) then
            do i = nelim, 1, -1
               ip = (i-1)*blkm
               rp1 = map(i)
               xlocal(i,r) = xlocal(i,r) - &
                  sum(xlocal(i+1:nelim,r) * lcol(ip+i+1:ip+nelim))
               xlocal(i,r) = xlocal(i,r) / lcol(ip+i)
               x(rp1,r) = xlocal(i,r)
            end do
         else
            do i = nelim, 1, -1
               ip = (i-1)*blkm
               rp1 = map(i)
               xlocal(i,r) = xlocal(i,r) - &
                  sum(xlocal(i+1:nelim,r) * lcol(ip+i+1:ip+nelim))
               x(rp1,r) = xlocal(i,r)
            end do
         end if
      end do
   end if
end subroutine slv_bwd_mult

!*************************************************************************
!
! Diagonal solve one rhs
!
subroutine slv_diag_one(invp, x, nelim, d, lperm)
   integer, dimension(*), intent(in) :: invp
   real(wp), dimension(*), intent(inout) :: x
   integer, intent(in) :: nelim
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm

   integer :: j
   integer :: rp1, rp2
   real(wp) :: temp

   j = 1
   do while(j .le. nelim)
      if (d(2*j).ne.0) then
         ! 2x2 pivot
         rp1 = invp( lperm(j) )
         rp2 = invp( lperm(j+1) )
         temp   = d(2*j-1) * x(rp1) + &
                  d(2*j)   * x(rp2)
         x(rp2) = d(2*j)   * x(rp1) + &
                  d(2*j+1) * x(rp2)
         x(rp1) = temp
         j = j + 2
      else
         ! 1x1 pivot
         rp1 = invp( lperm(j) )
         x(rp1) = x(rp1) * d(2*j-1)
         j = j + 1
      end if
   end do
end subroutine slv_diag_one

!*************************************************************************
!
! Diagonal solve multiple rhs
!
subroutine slv_diag_mult(invp, nrhs, x, ldx, nelim, d, lperm)
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,*), intent(inout) :: x
   integer, intent(in) :: nelim
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm

   integer :: j, r
   integer :: rp1, rp2
   real(wp) :: temp

   do r = 1, nrhs
      j = 1
      do while(j .le. nelim)
         if (d(2*j).ne.0) then
            ! 2x2 pivot
            rp1 = invp( lperm(j) )
            rp2 = invp( lperm(j+1) )
            temp     = d(2*j-1) * x(rp1,r) + &
                       d(2*j)   * x(rp2,r)
            x(rp2,r) = d(2*j)   * x(rp1,r) + &
                       d(2*j+1) * x(rp2,r)
            x(rp1,r) = temp
            j = j + 2
         else
            ! 1x1 pivot
            rp1 = invp( lperm(j) )
            x(rp1,r) = x(rp1,r) * d(2*j-1)
            j = j + 1
         end if
      end do
   end do
end subroutine slv_diag_mult

!*************************************************************************
!
! indefinite singular case. check consistency after forward solve
! and perform diag solve if consistent. 
!
subroutine check_diag(n, invp, nrhs, x, ldx, nelim, d, lperm, &
      flag_out, num_false, tol)
   integer, intent(in) :: n
   integer, dimension(*), intent(in) :: invp
   integer, intent(in) :: nrhs
   integer, intent(in) :: ldx
   real(wp), dimension(ldx,2*nrhs), intent(inout) :: x
   integer, intent(in) :: nelim
   real(wp), dimension(2*nelim) :: d
   integer, dimension(*), intent(in) :: lperm
   logical, intent (inout) :: flag_out(nrhs)
   integer, intent (inout) :: num_false ! number of inconsistent rhs
   real(wp), intent (in) :: tol(nrhs) ! tolerance used to determine whether
       ! or not system is consistent

   integer :: j, r
   integer :: rp1, rp2
   real(wp) :: temp

   ! write (6,*) 'nelim',nelim

   do r = 1, nrhs
      j = 1
      do while(j .le. nelim)

         if (d(2*j).ne.0) then
            ! 2x2 pivot
            rp1 = invp( lperm(j) )
            rp2 = invp( lperm(j+1) )
            temp     = d(2*j-1) * x(rp1,r) + &
                       d(2*j)   * x(rp2,r)
            x(rp2,r) = d(2*j)   * x(rp1,r) + &
                       d(2*j+1) * x(rp2,r)
            x(rp1,r) = temp
            j = j + 2
         else
            ! 1x1 pivot
            rp1 = invp( lperm(j) )
            temp = x(rp1,r)
            if (flag_out(r) .and. d(2*j-1).eq.zero) then
               ! check if rhs is inconsistent (but only if consistent so far)
               if (abs(temp) .gt. tol(r)) then
                  flag_out(r) = .false.
                  num_false = num_false + 1
                  ! set partial Fredholm solution to be zero,
                  ! except for entry corresponding to inconsistency
                  x(1:n,r+nrhs) = zero
                  x(rp1,r+nrhs) = one
               end if
            end if
            x(rp1,r) = temp * d(2*j-1)
            j = j + 1
         end if
      end do
   end do
end subroutine check_diag

!*************************************************************************
!
! Return diagonal entries to user
!
subroutine ma97_enquire_posdef_double(akeep, fkeep, control, info, d)
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), target, intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   real(wp), dimension(*), intent(out) :: d

   integer :: blkn, blkm
   character(50)  :: context      ! Procedure name (used when printing).
   integer(long) :: i
   integer :: j
   integer :: n
   integer :: node
   integer :: nout
   integer :: piv

   type(node_type), pointer :: nptr

   context = 'ma97_enquire_posdef' 
   info%flag = MA97_SUCCESS

   nout = control%unit_error
   if (control%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      ! immediate return if had an error
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (.not.fkeep%pos_def) then
      info%flag = MA97_ERROR_NOT_LLT
      call MA97_print_flag(context,nout,info%flag)
      return
   end if

   n = akeep%n
   ! ensure d is not returned undefined
   d(1:n) = zero ! ensure do not returned with this undefined
   
   piv = 1
   do node = 1, akeep%nnodes
      nptr => fkeep%nodes(node)
      blkn = akeep%sptr(node+1) - akeep%sptr(node)
      blkm = int(akeep%rptr(node+1) - akeep%rptr(node))
      i = 1
      do j = 1, blkn
         d(piv) = nptr%lcol(i)
         i = i + blkm + 1
         piv = piv + 1
      end do
   end do
end subroutine ma97_enquire_posdef_double

!*************************************************************************
! In indefinite case, the pivot sequence used will not necessarily be
!
! the same as that passed to ma97_factor (because of delayed pivots). This
! subroutine allows the user to obtain the pivot sequence that was
! actually used.
! also the entries of D^{-1} are returned using array d.
!
subroutine ma97_enquire_indef_double(akeep, fkeep, control, info, &
      piv_order, d)
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), target, intent(in) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info
   integer, dimension(*), optional, intent(out) :: piv_order
      ! If i is used to index a variable, its position in the pivot sequence
      ! will be placed in piv_order(i), with its sign negative if it is
      ! part of a 2 x 2 pivot; otherwise, piv_order(i) will be set to zero.
   real(wp), dimension(2,*), optional, intent(out) :: d ! The diagonal
      ! entries of D^{-1} will be placed in d(1,:i) and the off-diagonal
      ! entries will be placed in d(2,:). The entries are held in pivot order.

   integer :: blkn, blkm
   character(50)  :: context      ! Procedure name (used when printing).
   integer :: j, k
   integer :: n
   integer :: nd
   integer :: node
   integer :: nout
   integer(long) :: offset
   integer :: piv

   type(node_type), pointer :: nptr

   context = 'ma97_enquire_indef' 
   info%flag = MA97_SUCCESS

   nout = control%unit_error
   if (control%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      ! immediate return if had an error
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (fkeep%pos_def) then
      info%flag = MA97_ERROR_NOT_LDLT
      call MA97_print_flag(context,nout,info%flag)
      return
   end if

   n = akeep%n
   if(present(d)) then
      ! ensure d is not returned undefined
      d(1:2,1:n) = zero
   end if
   
   piv = 1
   do node = 1, akeep%nnodes
      nptr => fkeep%nodes(node)
      j = 1
      nd = nptr%ndelay
      blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
      blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
      offset = blkm*(blkn+0_long)
      do while(j .le. nptr%nelim)
         if (nptr%lcol(offset+2*j).ne.0) then
            ! 2x2 pivot
            if(present(piv_order))  then
               k = akeep%invp( nptr%perm(j) )
               piv_order(k) = -piv
               k = akeep%invp( nptr%perm(j+1) )
               piv_order(k) = -(piv+1)
            end if
            if(present(d)) then
               d(1,piv) = nptr%lcol(offset+2*j-1)
               d(2,piv) = nptr%lcol(offset+2*j)
               d(1,piv+1) = nptr%lcol(offset+2*j+1)
               d(2,piv+1) = 0
            end if
            piv = piv + 2
            j = j + 2
         else
            if(present(piv_order)) then
               k = akeep%invp( nptr%perm(j) )
               piv_order(k) = piv
            end if
            if(present(d)) then
               d(1,piv) = nptr%lcol(offset+2*j-1)
               d(2,piv) = 0
            end if
            piv = piv + 1
            j = j + 1
         end if
      end do
   end do
end subroutine ma97_enquire_indef_double

!*************************************************************************
!
! In indefinite case, the entries of D^{-1} may be changed using this routine.
!
subroutine ma97_alter_double(d, akeep, fkeep, control, info)
   real(wp), dimension(2,*), intent(in) :: d  ! The required diagonal entries
     ! of D^{-1} must be placed in d(1,i) (i = 1,...n)
     ! and the off-diagonal entries must be placed in d(2,i) (i = 1,...n-1).
   type(ma97_akeep), intent(in) :: akeep
   type(ma97_fkeep), target, intent(inout) :: fkeep
   type(ma97_control), intent(in) :: control
   type(ma97_info), intent(out) :: info

   integer :: blkm, blkn
   character(50)  :: context      ! Procedure name (used when printing).
   integer(long) :: ip
   integer :: j
   integer :: nd
   integer :: node
   integer :: nout
   integer :: piv

   type(node_type), pointer :: nptr

   info%flag = MA97_SUCCESS

   context = 'ma97_alter'

   nout = control%unit_error
   if (control%print_level < 0) nout = -1

   if (.not. allocated(fkeep%nodes)) then
      ! factorize phase has not been performed
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
    end if

   ! immediate return if already had an error
   if (akeep%flag .lt. 0 .or. fkeep%flag .lt. 0) then
      info%flag = MA97_ERROR_CALL_SEQUENCE
      call ma97_print_flag(context,nout,info%flag)
      return
   end if

   if (fkeep%pos_def) then
      info%flag = MA97_ERROR_NOT_LDLT
      call MA97_print_flag(context,nout,info%flag)
      return
   end if 
   
   piv = 1
   do node = 1, akeep%nnodes
      nptr => fkeep%nodes(node)
      nd = nptr%ndelay
      blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
      blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
      ip = blkm*(blkn+0_long) + 1
      do j = 1, nptr%nelim
         nptr%lcol(ip)   = d(1,piv)
         nptr%lcol(ip+1) = d(2,piv)
         ip = ip + 2
         piv = piv + 1
      end do
   end do
end subroutine ma97_alter_double

!*************************************************************************

subroutine free_akeep_double(akeep)
   type(ma97_akeep), intent(inout) :: akeep

   integer :: st

   deallocate(akeep%child_ptr, stat=st)
   deallocate(akeep%child_list, stat=st)
   deallocate(akeep%invp, stat=st)
   deallocate(akeep%level, stat=st)
   deallocate(akeep%nlist, stat=st)
   deallocate(akeep%nptr, stat=st)
   deallocate(akeep%rlist, stat=st)
   deallocate(akeep%rptr, stat=st)
   deallocate(akeep%sparent, stat=st)
   deallocate(akeep%sptr, stat=st)
   deallocate(akeep%subtree_work, stat=st)
   deallocate(akeep%ptr, stat=st)
   deallocate(akeep%row, stat=st)
   deallocate(akeep%map, stat=st)
end subroutine free_akeep_double

!*******************************

subroutine free_fkeep_double(fkeep)
   type(ma97_fkeep), intent(inout) :: fkeep

   integer :: st

   if (.not.allocated(fkeep%nodes)) return

   call smfreeall(fkeep%alloc)
   deallocate(fkeep%alloc)
   nullify(fkeep%alloc)

   deallocate(fkeep%nodes, stat=st)
   deallocate(fkeep%scaling, stat=st)
end subroutine free_fkeep_double

!*******************************

subroutine finalise_both_double(akeep, fkeep)
   type(ma97_akeep), intent(inout) :: akeep
   type(ma97_fkeep), intent(inout) :: fkeep

   call free_akeep_double(akeep)
   call free_fkeep_double(fkeep)
end subroutine finalise_both_double

!***************************************************************
!
! This subroutine recursivly subdivides factorization at a node
! ( L_11      )
! ( L_21 L_22 )
! 1. Any delays are swapped to the end
! 2. Column L_:1 is factorized (recusive call)
! 3. Update L_22 using L_21
! 4. Swap any delays from L_:1 to end
! 5. Factorize L_22 (recursive call)
! 6. Update any delays from L_:1 by L_22.
! 7. Try an eliminate any delays again.
!
! If the width of the column is smaller than minsz then the factorization
! kernel factor_solve_block is called.
recursive subroutine rfact_block(pos_def, n, p, nb, s, a, lda, d, buf, &
      perm,  q, control, stats, aleft, aln, ntask, sublevel)
   logical, intent(in) :: pos_def
   integer, intent(in) :: n ! number of rows in trapeziodal matrix
   integer, intent(in) :: p ! number of cols in trapeziodal matrix
   integer, intent(in) :: nb ! Block size
   integer, intent(in) :: s ! number of cols in trapeziodal matrix that
      ! are from delayed pivots. These cols are searched last for pivots.
   integer :: lda  ! leading dimension of a
   real(wp), intent(inout) :: a(*) ! holds trapezoidal matrix 
      ! to be factorized. Holds the factorization on exit.
      ! Delayed columns are permuted to the final columns. Held rowwise
      ! at all times.
   real(wp), intent(inout) :: d(*) ! holds diagonal pivots
   type(real_ptr_type), dimension(*), intent(inout) :: buf ! Work array
   integer, intent(inout) :: perm(*) ! perm array. Any swaps performed
      ! on columns of matrix are also performed on this array.
   integer, intent(out) :: q ! number of pivots chosen
   type(ma97_control), intent(in) :: control
   type(thread_stats), dimension(*), intent(inout) :: stats
   real(wp), intent(inout) :: aleft(*)
   integer, intent(in) :: aln
   integer, intent(in) :: ntask
   integer, optional, intent(in) :: sublevel

   integer :: ndelay
   integer :: i, j
   integer :: mu, nu
   integer :: en2
   integer :: q1
   integer(long) :: ip, ip2
   integer :: lapack_info
   integer :: this_thread
   integer :: maxdepth

   integer, parameter :: minsz = 256 ! minimum size to engage recursion
   !%%type(log_type) :: ltask

   if(present(sublevel)) then
      maxdepth = sublevel-1
   elseif(p.le.minsz) then
      maxdepth = 1 ! Avoid expensive log operations
   else
      maxdepth = max(1, int(log(real(n/minsz)) / log(2.0)))
   end if

   this_thread = 1
!$ this_thread = omp_get_thread_num() + 1

   if (p .le. minsz .or. maxdepth.eq.1) then
      if(pos_def) then
         call dpotrf('L', p, a, lda, lapack_info)
         if(lapack_info.ne.0) then
            stats(this_thread)%flag = MA97_ERROR_NOT_POS_DEF
            return
         end if
         call dtrsm('Right', 'Lower', 'Trans', 'Non-unit', n-p, p, &
            one, a, lda, a(p+1), lda)
         q = p
      else
         nullify(buf(this_thread)%chkptr)
         if(size(buf(this_thread)%val) .lt. p**2) then
            deallocate(buf(this_thread)%val, stat=stats(this_thread)%st)
            allocate(buf(this_thread)%val(p**2), stat=stats(this_thread)%st)
            if(stats(this_thread)%st.ne.0) then
               stats(this_thread)%flag = MA97_ERROR_ALLOCATION
               return
            end if
         end if
         call factor_solve_block(n, p, nb, s, a, lda, d, &
            buf(this_thread)%val, perm, q, control, stats(this_thread), &
            aleft, aln, recheck=.not.present(sublevel))
      end if
   else
      ! Swap any delays we have been passed to end
      if (s>0 .and. s<p) then
         i = min(s,p-s)
         ! Make first i columns be last
         do j = 1,i
            nullify(buf(this_thread)%chkptr)
            call ma64_swap(n,nb,0,1,0,a,lda,aleft,aln, &
               buf(this_thread)%val,p,perm,j,p+1-j)
         end do
      end if

      ! Factor first half
      en2 = p/2
      call rfact_block(pos_def, n, en2, nb, 0, a, lda, d, buf, perm, q1, &
         control, stats, aleft, aln, ntask, sublevel=maxdepth)
      if(stats(this_thread)%flag.lt.0) return

      ! Update second half
      ! Note: all first half columns are already up to date, but only first q
      ! should be used for update!
      ip = en2+1
      ip2 = en2*lda + en2+1
      mu = n - en2
      nu = p - en2
      if (q1.ge.1) then
         if(pos_def) then
            call mydsyrk(mu, nu, q1, -one, a(ip), lda, &
               one, a(ip2), lda, control%min_ldsrk_work)
         else
            call ldsrk(mu, nu, q1, -one, a(ip), lda, d, &
               one, a(ip2), lda, buf, control%min_ldsrk_work, stats)
            if(stats(this_thread)%flag.lt.0) return
         end if
      end if

      ! Push delays to end
      ! Note: this means that some previously rejected pivots from child
      ! nodes may get moved into the middle again. They will hopefully be better
      ! pivot candidates than those they are replaced by however.
      ndelay = en2 - q1
      if (ndelay.gt.0) then
         ! Swap column set [q1+1:en2] to end
         do j = q1+1, en2
            nullify(buf(this_thread)%chkptr)
            call ma64_swap(n,nb,0,1,0,a,lda,aleft,aln, &
               buf(this_thread)%val,1,perm,j,p-en2+j)
         end do
      end if

      ! Factor second half (includes delays from first half)
      ip = q1*lda + q1 + 1
      call rfact_block(pos_def, n-q1, p-q1, nb, 0, a(ip), lda, &
         d(2*q1+1), buf, perm(q1+1), q, control, stats, aleft(q1+1), &
         aln+q1, ntask, sublevel=maxdepth)
      if(stats(this_thread)%flag.lt.0) return

      ! Update pivoted upon count
      q1 = q + q1

      ! If we are the top level, force through any delays we can by iterating
      ! over delays repeatedly until no more can be eliminated.
      q = 0
      if (p-q1.gt.0 .and. .not.present(sublevel)) then
         ip = q1*lda + q1 + 1
         nullify(buf(this_thread)%chkptr)
         if(size(buf(this_thread)%val) .lt. (p-q1)**2) then
            deallocate(buf(this_thread)%val, stat=stats(this_thread)%st)
            allocate(buf(this_thread)%val((p-q1)**2), &
               stat=stats(this_thread)%st)
            if(stats(this_thread)%st.ne.0) then
               stats(this_thread)%flag = MA97_ERROR_ALLOCATION
               return
            end if
         end if
         call factor_solve_block(n-q1, p-q1, nb, 0, a(ip), lda, &
            d(2*q1+1), buf(this_thread)%val, perm(q1+1), q, control, &
            stats(this_thread), aleft(q1+1), aln+q1, recheck=.true.)
      end if

      q = q + q1
   end if

end subroutine rfact_block

!*******************************
!
! performs factorization of trapezoidal matrix.
! based on hsl_ma64 version 6.0.0 (ported 18th January 2011)
!
! JAS: 23 January 2012. Minor tidying of code (removal of redundant
! variables, added comments etc)
!
subroutine factor_solve_block(n, p, nb, s, a, lda, d, buf, perm,  &
      q, control, stats, aleft, aln, recheck)
   integer, intent(in) :: n ! number of rows in trapeziodal matrix
   integer, intent(in) :: p ! number of candidate pivots
   integer, intent (in) :: nb ! Column block size
   integer, intent(in) :: s ! number of cols in trapeziodal matrix that
      ! are from delayed pivots. These cols are swapped to the
      ! end and are searched last for pivots.
   integer, intent(in) :: lda  ! leading dimension of a
   real(wp), intent(inout) :: a(*) ! holds trapezoidal matrix 
      ! to be factorized. Holds the factorization on exit.
      ! Delayed columns are permuted to the final columns.
   real(wp), intent(inout) :: d(2*p)
   real (wp) :: buf(p*p) ! Work array
   integer, intent(inout) :: perm(*) ! perm array. Any swaps performed
      ! on columns of matrix are also performed on this array.
   integer, intent(out) :: q ! number of pivots chosen
   type(ma97_control), intent(in) :: control
   type(thread_stats), intent(inout) :: stats
   real(wp) :: aleft(*) ! set of rows that need swapping too, leading dim lda
   integer :: aln ! number of columns in aleft
   logical, intent(in) :: recheck ! controls whether to check through the
      ! columns more than once while not all columns pivoted on

   !     .. Local Scalars ..
   real(wp) abs_amax ! Largest absolute value to left of diagonal in row m.
   real(wp) amax ! Entry of largest absolute value to left of diagonal in
      ! row m.
   real(wp) amax2 ! Second largest absolute value to left of diagonal in row m
   real(wp) amaxm1 ! Largest absolute value below diagonal in column m, not
      ! including row m+1
   real(wp) amaxb ! Largest absolute value below diagonal in column m
   real(wp) amaxt ! Largest absolute value in column t.
   real(wp) amaxt_cache ! Largest absolute value in column m-1 (used as cache).
   real(wp) aval ! Temporary work variable
   real(wp) detpiv ! Determinant of candidate pivot is detpiv/detscale.
   real(wp) detpiv0 ! First term in calculation of detpiv.
   real(wp) detpiv1 ! Second term in calculation of detpiv.
   real(wp) detscale ! Inverse of the largest entry in the candidate pivot.
   integer i ! Row index
   integer(long) ii ! long integer loop variable
   integer j ! Column index
   integer(long) jj ! long integer loop variable
   integer(long) k  ! Position in a
   integer(long) kb1 ! Position in buf of diagonal of column q+1
   integer(long) kkj ! Position of diagonal of column j
   integer(long) kkt ! Position of diagonal of column t
   integer(long) kkr1 ! Position of diagonal of column r+1
   integer(long) kkq !  Start of the pivotal block column
   integer(long) kkq1 ! Position of diagonal of column q+1
   integer(long) kkm ! Position of diagonal of column m
   integer(long) kq1 ! Position in a of diagonal of column q+1
   integer ldbuf ! leading dimension of buffer. It is set to p
   integer m  ! Column searched most recently
   integer mbest ! Column with best relative pivot value for a 1x1 pivot
   integer mdummy  ! Loop execution count
   integer mlast ! Last column of the block in which m appears
   integer nbi ! Inner block size. This is set to 16.
   integer pivsiz ! Size of the chosen pivot, 0 if none chosen, or -1
      ! if column is essentially zero
   real(wp) pivval ! temporary variable for storing pivotal value
   integer qlast ! Last column of the inner block containing column q+1
   integer rm ! Number of pivot operations applied to the (outer) block
              ! containing column m
   real(wp) rmax ! Largest entry in column m
   real(wp) rmax2 ! Largest entry in column m outwith rows m,m-1.
   integer t ! Candidate 2x2 pivot is in columns t and m
   real(wp) u ! Relative pivot threshold
   real(wp) ubest1 ! Relative pivot value of best candidate 1x1 pivot
   real(wp) urel ! Relative pivot value
   integer :: nzero ! Number of zero pivots since last non-zero pivot

   logical :: noreset

   ! Treat special (small) cases first.
   q = 0
   if (p.eq.0) return

   if (p.eq.1) then
      call factor_solve_block1(n, a, d, q, control, stats)
      return
   end if

   if(n*p < 3072) then ! fits within 24KB
      call factor_solve_block_small(n, p, nb, s, a, lda, d, buf, perm,  &
         q, control, stats, aleft, aln, recheck)
      return
   end if

   u = min(max(control%u,zero),one)

   nbi = 16

   kkq1 = 1
   kkq = 1
   m = p
   qlast = min(nbi,p)
   ldbuf = p ! not changed during calculation, just used to make code more
             ! maintainable
   ! m is updated at the start of the main loop so initializing it to p causes
   ! it to have the value 1 during the first execution of that loop.

   !print *
   !print *, "============================================="
   !print *, "factor n, p, s = ", n, p, s
   !print *, "a input = "
   !print "(11es12.4)", a

   if (s>0 .and. s<p) then
      i = min(s,p-s)
      ! Make first i columns be last
      do j = 1,i
         call ma64_swap(n,nb,q,1,0,a,lda,aleft,aln,buf,ldbuf,perm,j,p+1-j)
      end do
   end if

   nzero = 0 ! number of zero columns encountered since last proper pivot
   noreset = .false.
   pivot: do
      !print *, ".eq. q=", q, ".eq."
      ! Perform a pivotal operation
      pivsiz = 0
      ubest1 = zero
      amaxt_cache = -1 ! cache is intially invalid
      sweep: do mdummy = 1, p-q ! Look for a pivotal column or pair of columns
         ! Update m and the scalars associated with column m.
         ! Remember: on first time through, m = p
         m = m+1
         if (m>p) then
            if(noreset) exit sweep
            ! Go back to column q+1
            m = q+1
            kkm = kkq1
            rm = q
            noreset = .not.recheck
         else if (m<1+nb) then
            ! Within the current block column
            kkm = kkm + lda + 1
         end if

         ! Update column m by pivots so far selected within inner block
         kkr1 = kkq1 + (rm-q)*(lda+1_long)
         k = kkr1+m-rm-1
         if (q>rm) then
            !call dgemv('NoTrans',n-m+1,q-rm,-one,a(k),lda, &
            !   buf(lda*(rm+0_long)+m),lda,one,a(kkm),1)
            call mydgemv(n-m+1,q-rm,a(k),lda, &
               buf(ldbuf*(rm+0_long)+m),ldbuf,a(kkm))
         end if

         !print *, "scan ", m
         ! Find largest entry (abs_amax) and second largest entry (amax2)
         ! to left of diagonal in row m (cols q+1:m-1).
         ! Let t be row index for abs_amax
         j = q + 1
         k = kkq1 + m - j
         amax = zero
         abs_amax = zero
         amax2 = zero
         t = 0
         if (j<m) then
            t = j
            kkt = kkq1
         end if
         do j = j, m-1
            if (abs(a(k))>abs_amax) then
               t = j
               amax2 = abs(amax)
               amax = a(k)
               abs_amax = abs(a(k))
               kkt = k - (m-j)
            else
               !amax2 = max(abs(a(k)),amax2)
               if (abs(a(k)).gt.amax2) amax2 = abs(a(k))
            end if
            k = k + lda
         end do

         ! Now calculate largest entry below the diagonal of column m.
         amaxm1 = zero
         amaxb = zero
         !do i = m+1,n
         !   amaxb = max(abs(a(kkm+i-m)),amaxb)
         !end do
         if (n-m.ge.1) amaxb = abs(a(kkm+1))
         do k = kkm+2,kkm+n-m
            if (abs(a(k)).gt.amaxm1) amaxm1 = abs(a(k))
         end do
         if (amaxm1.gt.amaxb) amaxb = amaxm1

         ! Now calculate largest entry in the whole of column m and make sure
         ! that it is neither small nor infinity.
         !rmax = max(abs_amax,abs(a(kkm)),amaxb)
         rmax = abs_amax
         if (abs(a(kkm)).gt.rmax) rmax = abs(a(kkm))
         if (amaxb.gt.rmax) rmax = amaxb
         if (rmax<=control%small) then
            ! All entries of the column are small
            a(kkm) = zero
            pivsiz = -1
            perm(m) = perm(m)
            exit sweep
         else if (rmax > huge(zero)) then
            ! There is an infinity in the column
            stats%flag = MA97_ERROR_INFINITY
            return
         end if

         ! Calculate the relative pivot value and see if it is the best so far
         if (abs(a(kkm))>control%small) then
            urel = abs(a(kkm))/rmax
         else
            urel = zero
         end if
         if (urel >= ubest1) then
            ubest1 = urel
            mbest = m
         end if

         ! If there is a candidate 2x2 pivot, try it.
         ! Look for a 1x1 pivot only if the 2x2 pivot is unacceptable.
         tgt0: if (t>0) then
            if ( min(abs(a(kkm)),abs(a(kkt))).gt.control%small .or. &
                  abs_amax.ge.control%small) then

               ! Store value of the largest entry in whole of column m outwith
               ! rows m and t
               rmax2 = max(amax2,amaxb)

               detscale = one/max( abs(a(kkm)), abs(a(kkt)), abs_amax )
               detpiv1 =  (amax*detscale)*amax
               detpiv0 = a(kkm)*detscale*a(kkt)
               detpiv = detpiv0 - detpiv1
               ! Make sure that the 2x2 pivot is not singular and that there is
               ! little cancellation in calculating detpiv. Bearing in mind the
               ! way detscale is calculated, if the largest entry of the matrix
               ! is chosen as pivot, the one entry of the reduced matrix has
               ! absolute value abs(detpiv).
               !print *, "test 2x2 = ", detpiv, detpiv0, detpiv1
               left2x2:if (abs(detpiv)> &
                     max(control%small,abs(detpiv0)/2,abs(detpiv1)/2)) then

                  ! Find largest entry in column t outwith rows m and t
                  if (t.eq.m-1 .and. amaxt_cache.ne.-1) then
                     ! Use cached answer from scan of previous column
                     amaxt = amaxt_cache
                  else
                     amaxt = zero
                     j = q + 1
                     k = kkq1 + t - j - lda
                     do j = q+1, t-1
                        k = k + lda
                        amaxt = max(abs(a(k)),amaxt)
                     end do
                     k = k + lda
                     do i = t+1,m-1
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                     do i = m+1,n
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                  end if

                  ! OK as 2x2 pivot if all entries in the rest of the columns
                  ! are small
                  if (max(rmax2,amaxt)<=control%small) then
                     pivsiz = 2
                     exit sweep
                  end if

                  ! Calculate the relative pivot value (as 2x2 pivot)
                  urel = abs(detpiv)/max( &
                     abs(a(kkm)*detscale)*amaxt+(abs_amax*detscale)*rmax2, &
                     abs(a(kkt)*detscale)*rmax2+(abs_amax*detscale)*amaxt )
                  !print *, "urel2x2 = ", urel

                  ! OK as 2x2 pivot if relative pivot value is big enough
                  if (urel>u) then
                     pivsiz = 2
                     exit sweep
                  end if
               end if left2x2
            end if
         end if tgt0
         amaxt_cache = max(abs_amax, amaxm1)

         ! If 2x2 pivot rejected or only one column left, take best 1x1 pivot
         ! if it is OK.
         if (t>0 .or. m.eq.p) then
            !print *, "   test 1x1 ubest1 = ", ubest1
            if (ubest1>u) then
               !if (t>0) print *, "accept 1x1 rather than 2x2"
               pivsiz = 1
               if (mbest.ne.m) &
                  call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf, &
                     perm,mbest,m)
               exit sweep
            end if
         end if
      end do sweep

      !
      ! At this stage following variables should have good values:
      ! m - pivot column (wil be swapped to posn q+1 for 1x1 or q+2 for 2x2)
      ! q - number of pivots already performed
      ! kkq1 - posn in a of diagonal of col q+1
      ! t - other pivot column for 2x2 (t<m) (swapped to posn q+q
      ! detscale & detpiv - used for stats
      !

      pivsiz0: if (pivsiz.eq.0) then
         !print *, "pivsize0 dropout", ubest1
         ! No pivot found in search of all available columns
         ! Since all the columns have been updated, revise m to q+1
         m = q+1
         kkm = kkq1
         rm = q
      end if pivsiz0

      pivsizes: if (pivsiz.eq.1) then
         !print *, "1x1 pivot on ", m
         nzero = 0 ! pivot found, reset recent zero count
         ! Swap columns if m not q+1
         if (q+1.ne.m) &
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm,q+1,m)
         ! We store D**-1.
         !print *, "1x1 pivot ", a(kkq1), "(col ", m, ")"
         d(2*q+1) = one/a(kkq1)
         d(2*q+2) = zero
         if (a(kkq1).lt.0) stats%num_neg = stats%num_neg + 1
         ! Store L in A and LD in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = a(kkq1)
         a(kkq1) = 1/a(kkq1)
         pivval = d(2*q+1)
         jj = kb1+1
         do ii = kq1+1, kq1+p-q-1
            buf(jj) = a(ii); jj = jj + 1
            a(ii) = pivval*a(ii)
         end do
         do ii = kq1+p-q, kq1+n-q-1
            a(ii) = pivval*a(ii)
         end do
         ! Update columns q+2 to m (ie keep up to date columns that have 
         ! already been tested but failed)
         kkj = kkq1
         do j = q+2, min(nb,m)
            kkj = kkj + lda + 1
            call daxpy(n-j+1,-buf(kb1+j-q-1),a(kkq1+j-q-1),1,a(kkj),1)
         end do
         ! Update q and kkq1
         kkq1 = kkq1 + lda + 1
         q = q+1

      else if (pivsiz.eq.2) then pivsizes
         !print *, "2x2 pivot on ", t, m
         nzero = 0 ! pivot found, reset recent zero count
         ! Swap columns unless t.eq.q+1 and m.eq.q+2
         if (q+2.ne.m) then
            if (q+1.ne.t) &
               call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf, &
                  perm,q+1,t)
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm, &
               q+2,m)
         end if
         ! We store D**-1
         k = kkq1 + lda + 1
         d(2*q+1) = (a(k)*detscale)/detpiv
         d(2*q+3) = (a(kkq1)*detscale)/detpiv
         d(2*q+2) = -(a(kkq1+1)*detscale)/detpiv
         d(2*q+4) = zero
         !print *, "2x2 pivot ", d(2*q+1:2*q+3)
         !print *, "   (cols ", t, m, ")"
         !print *, "   [vals ", a(k), a(kkq1), a(kkq1+1), &
         !   "dp", detscale, detpiv, "]"
         ! Update info
         if (detpiv<zero) then
            stats%num_neg = stats%num_neg + 1
         else if (a(kkq1)+a(k)<0) then
            stats%num_neg = stats%num_neg + 2
         end if
         stats%num_two = stats%num_two + 1
         ! Store L in A and LD or conjg(LD) in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = a(kkq1)
         buf(kb1+1) = a(kkq1+1)
         buf(kb1+ldbuf) = a(kkq1+1)
         buf(kb1+ldbuf+1) = a(k)
         a(kkq1) = 1.0_wp
         a(k) = 1.0_wp
         a(kkq1+1) = 0.0_wp
         do i = 2, p-q-1
            buf(kb1+i      ) = a(kkq1+i)
            buf(kb1+i+ldbuf) = a(k+i-1)
            a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
            a(k+i-1) = d(2*q+3)*a(k+i-1) + d(2*q+2)*buf(kb1+i)
         end do
         do i = p-q, n-q-1
            aval = a(kkq1+i)
            a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
            a(k+i-1) = d(2*q+3)*a(k+i-1) + d(2*q+2)*aval
         end do
         ! Update columns q+3 to m (ie keep up to date columns that have 
         ! already been tested but failed)
         kkj = k + lda + 1
         j = max(q+3,1)
         if (m>=j) then
            call dgemm('n','t',n-j+1,m-j+1,2,-one,a(kkq1+j-q-1), &
               lda-((q+1)/nb)*nb,buf(kb1+j-q-1),ldbuf,one,a(kkj),lda)
         end if
         ! Update q and kkq1
         kkq1 = k + lda + 1
         q = q + 2

      else if (pivsiz.eq.-1) then pivsizes
         !print *, "ZERO pivot"
         nzero = nzero + 1 ! count zero pivot
         ! Handle a row that is zero
         if (q+1.ne.m) &
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm, &
               q+1,m)
         !print *, "fail piv"
         d(2*q+1) = zero
         d(2*q+2) = zero
         if(.not.control%action) then
            ! Error on singular matrix
            stats%flag = MA97_ERROR_SINGULAR
            return
         end if
         stats%num_zero = stats%num_zero + 1
         ! Store L in A and -LD in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = zero
         a(kkq1) = zero
         do i = 1, p-q-1
            buf(kb1+i) = zero
            a(kkq1+i) = zero
         end do
         do i = p-q, n-q-1
            a(kkq1+i) = zero
         end do
         kkq1 = kkq1 + lda + 1
         q = q+1
      end if pivsizes

      if (q>=qlast) then
         ! Inner block is complete
         ! Update columns (m+1:min(p,m1+nb-1)) for pivots rm+1:qlast (BLAS3)
         mlast = min(p,nb)
         if (m<mlast) call ma64_update (n,nb,kkq,a,lda,buf,ldbuf,m+1, &
            mlast,rm+1,qlast-nzero)
         rm = qlast
         qlast = min (p,qlast+nbi)
         nzero = 0 ! reset recent zero count
      end if

      if (q .eq. p .or. pivsiz.eq.0) exit pivot
   end do pivot

   !print *, "exit d = ", d(1:2*q)
   !print *, "exit a = "
   !print "(5es12.4)", a
   !print *, "delayed ", p-q, "pivots"

end subroutine factor_solve_block

!*******************************
!
! This version when full block column easily fits within L1 cache
! We don't use BLAS, and implement a simple right-looking kernel
!
subroutine factor_solve_block_small(n, p, nb, s, a, lda, d, buf, perm,  &
      q, control, stats, aleft, aln, recheck)
   integer, intent(in) :: n ! number of rows in trapeziodal matrix
   integer, intent(in) :: p ! number of cols in trapeziodal matrix
   integer, intent (in) :: nb ! Block size
   integer, intent(in) :: s ! number of cols in trapeziodal matrix that
     ! are from delayed pivots. These cols are searched last for pivots.
   integer, intent(in) :: lda  ! leading dimension of a
   real(wp), intent(inout) :: a(*) ! holds trapezoidal matrix 
     ! to be factorized. Holds the factorization on exit.
      ! Delayed columns are permuted to the final columns. Held rowwise
      ! at all times.
   real(wp), intent(inout) :: d(2*p)
   real (wp) :: buf(p*p) ! Work array
   integer, intent(inout) :: perm(*) ! perm array. Any swaps performed
      ! on columns of matrix are also performed on this array.
   integer, intent(out) :: q ! number of pivots chosen
   type(ma97_control), intent(in) :: control
   type(thread_stats), intent(inout) :: stats
   real(wp) :: aleft(*) ! set of rows that need swapping too, leading dim lda
   integer :: aln ! number of columns in aleft
   logical, intent(in) :: recheck

   !     .. Local Scalars ..
   real(wp) abs_amax ! Largest absolute value to left of diagonal in row m.
   real(wp) amax ! Entry of largest absolute value to left of diagonal in
      ! row m.
   real(wp) amax2 ! Second largest absolute value to left of diagonal in row m
   real(wp) amaxm1 ! Largest absolute value below diagonal in column m, not
      ! including row m+1
   real(wp) amaxb ! Largest absolute value below diagonal in column m
   real(wp) amaxt ! Largest absolute value in column t.
   real(wp) amaxt_cache ! Largest absolute value in column m-1 (used as cache).
   real(wp) aval ! Temporary work variable
   real(wp) detpiv ! Determinant of candidate pivot is detpiv/detscale.
   real(wp) detpiv0 ! First term in calculation of detpiv.
   real(wp) detpiv1 ! Second term in calculation of detpiv.
   real(wp) detscale ! Inverse of the largest entry in the candidate pivot.
   integer i ! Row index
   integer(long) ii ! long integer loop variable
   integer j ! Column index
   integer(long) jj ! long integer loop variable
   integer(long) k  ! Position in a
   integer(long) kb1 ! Position in buf of diagonal of column q+1
   integer(long) kkj ! Position of diagonal of column j
   integer(long) kkt ! Position of diagonal of column t
   integer(long) kkq1 ! Position of diagonal of column q+1
   integer(long) kkm ! Position of diagonal of column m
   integer(long) kq1 ! Position in a of diagonal of column q+1
   integer ldbuf ! leading dimension of buffer. Should be set ldbuf=p
   integer m  ! Column searched most recently
   integer mbest ! Column with best relative pivot value for a 1x1 pivot
   integer mdummy  ! Loop execution count
   integer pivsiz ! Size of the chosen pivot, 0 if none chosen, or -1
      ! if column is essentially zero
   real(wp) pivval ! temporary variable for storing pivotal value
   integer rm ! Number of pivot operations applied to the (outer) block
              ! containing column m
   real(wp) rmax ! Largest entry in column m
   real(wp) rmax2 ! Largest entry in column m outwith rows m,m-1.
   integer t ! Candidate 2x2 pivot is in columns t and m
   real(wp) u ! Relative pivot threshold
   real(wp) ubest1 ! Relative pivot value of best candidate 1x1 pivot
   real(wp) urel ! Relative pivot value
   real(wp) umin ! Minimum relative pivot threshold
   integer :: nzero ! Number of zero pivots since last non-zero pivot

   real(wp) :: alpha
   real(wp) :: d11, d22

   logical :: noreset

   ! u is reset to control%u for each block column, but may then be relaxed to
   ! umin on a block column by block column basis
   u = min(max(control%u,zero),one)

   umin = u
   if (p.eq.n) umin = min(umin,0.5_wp)
   q = 0
   if (p.eq.0) return
   kkq1 = 1
   m = p
   ldbuf = p ! not changed during calculation, just used to make code more
      ! maintainable
   ! m is updated at the start of the main loop so initializing it to p causes
   ! it to have the value 1 during the first execution of that loop.


   !print *
   !print *, "============================================="
   !print *, "factor n, p, s = ", n, p, s
   !print *, "a input = "
   !print "(11es12.4)", a

   if (s>0 .and. s<p) then
      i = min(s,p-s)
      ! Make first i columns be last
      do j = 1,i
         call ma64_swap(n,nb,q,1,0,a,lda,aleft,aln,buf,ldbuf,perm,j,p+1-j)
      end do
   end if

   nzero = 0 ! number of zero columns encountered since last proper pivot
   noreset = .false.
   pivot: do
      !print *, ".eq. q=", q, ".eq."
      ! Perform a pivotal operation
      pivsiz = 0
      ubest1 = zero
      amaxt_cache = -1 ! cache is intially invalid
      sweep: do mdummy = 1, p-q ! Look for a pivotal column or pair of columns
         ! Update m and the scalars associated with column m
         m = m+1
         if (m>p) then
            if(noreset) exit sweep
            ! Go back to column q+1
            m = q+1
            kkm = kkq1
            rm = q
            noreset = .not.recheck
         else if (m<1+nb) then
            ! Within the current block column
            kkm = kkm + lda + 1
         end if

         !print *, "scan ", m
         ! Find largest and second largest entry to left of diagonal in row m.
         j = q + 1
         k = kkq1 + m - j - lda
         amax = zero
         abs_amax = zero
         amax2 = zero
         t = 0
         if (j<m) then
            t = j
            kkt = kkq1
         end if
         do j = j, m-1
            k = k + lda
            if (abs(a(k))>abs_amax) then
               t = j
               amax2 = abs(amax)
               amax = a(k)
               abs_amax = abs(a(k))
               kkt = k - (m-j)
            else
               !amax2 = max(abs(a(k)),amax2)
               if (abs(a(k)).gt.amax2) amax2 = abs(a(k))
            end if
         end do
         rmax = abs_amax
         if (abs(a(kkm)).gt.rmax) rmax = abs(a(kkm))

         ! Now calculate largest entry below the diagonal of column m.
         amaxm1 = zero
         amaxb = zero
         !do i = m+1,n
         !   amaxb = max(abs(a(kkm+i-m)),amaxb)
         !end do
         if (n-m.ge.1) amaxb = abs(a(kkm+1))
         do k = kkm+2,kkm+n-m
            if (abs(a(k)).gt.amaxm1) amaxm1 = abs(a(k))
         end do
         if (amaxm1.gt.amaxb) amaxb = amaxm1

         ! Now calculate largest entry in the whole of column m and make sure
         ! that it is neither small nor infinity.
         !rmax = max(abs_amax,abs(a(kkm)),amaxb)
         if (amaxb.gt.rmax) rmax = amaxb
         if (rmax<=control%small) then
            ! All entries of the column are small
            a(kkm) = zero
            pivsiz = -1
            perm(m) = perm(m)
            exit sweep
         else if (rmax > huge(zero)) then
            ! There is an infinity in the column
            stats%flag = MA97_ERROR_INFINITY
            return
         end if

         d11 = a(kkm)

         ! Calculate relative pivot value and see if it is the best so far
         if (abs(d11)>control%small) then
            urel = abs(d11)/rmax
            if(urel.gt.u) then
               pivsiz = 1
               exit sweep
            end if
            !if (urel >= ubest1) then
            !   ubest1 = urel
            !   mbest = m
            !end if
         end if

         ! If there is a candidate 2x2 pivot, try it.
         ! Look for a 1x1 pivot only if the 2x2 pivot is unacceptable.
         tgt0: if (t>0) then
            d22 = a(kkt)
            if ( min(abs(d11),abs(d22)).gt.control%small .or. &
                  abs_amax.ge.control%small) then

               ! Store value of the largest entry in whole of column m outwith
               ! rows m and t
               rmax2 = max(amax2,amaxb)

               detscale = one/max( abs(d11), abs(d22), abs_amax )
               detpiv1 =  (amax*detscale)*amax
               detpiv0 = d11*detscale*d22
               detpiv = detpiv0 - detpiv1
               ! Make sure that the 2x2 pivot is not singular and that there is
               ! little cancellation in calculating detpiv. Bearing in mind the
               ! way detscale is calculated, if the largest entry of the matrix
               ! is chosen as pivot, the one entry of the reduced matrix has
               ! absolute value abs(detpiv).
               !print *, "test 2x2 = ", detpiv, detpiv0, detpiv1
               left2x2:if (abs(detpiv)> &
                     max(control%small,abs(detpiv0)/2,abs(detpiv1)/2)) then

                  ! Find largest entry in column t outwith rows m and t
                  if (t.eq.m-1 .and. amaxt_cache.ne.-1) then
                     ! Use cached answer from scan of previous column
                     amaxt = amaxt_cache
                  else
                     amaxt = zero
                     j = q + 1
                     k = kkq1 + t - j - lda
                     do j = q+1, t-1
                        k = k + lda
                        amaxt = max(abs(a(k)),amaxt)
                     end do
                     k = k + lda
                     do i = t+1,m-1
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                     do i = m+1,n
                        amaxt = max(abs(a(k+i-t)),amaxt)
                     end do
                  end if

                  ! OK as 2x2 pivot if all entries in the rest of the columns
                  ! are small
                  if (max(rmax2,amaxt)<=control%small) then
                     pivsiz = 2
                     exit sweep
                  end if

                  ! Calculate the relative pivot value (as 2x2 pivot)
                  urel = abs(detpiv)/max( &
                     abs(d11*detscale)*amaxt+(abs_amax*detscale)*rmax2, &
                     abs(d22*detscale)*rmax2+(abs_amax*detscale)*amaxt )
                  !print *, "urel2x2 = ", urel

                  ! OK as 2x2 pivot if relative pivot value is big enough
                  if (urel>u) then
                     pivsiz = 2
                     exit sweep
                  end if
               end if left2x2
            end if
         end if tgt0
         amaxt_cache = max(abs_amax, amaxm1)

         ! If 2x2 pivot rejected or only one column left, take best 1x1 pivot
         ! if it is OK.
         if (t>0 .or. m.eq.p) then
            !print *, "   test 1x1 ubest1 = ", ubest1
            if (ubest1>u) then
               !if (t>0) print *, "accept 1x1 rather than 2x2"
               pivsiz = 1
               if (mbest.ne.m) &
                  call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm, &
                     mbest,m)
               exit sweep
            end if
         end if
      end do sweep

      !
      ! At this stage following variables should have good values:
      ! m - pivot column (wil be swapped to posn q+1 for 1x1 or q+2 for 2x2)
      ! q - number of pivots already performed
      ! kkq1 - posn in a of diagonal of col q+1
      ! t - other pivot column for 2x2 (t<m) (swapped to posn q+q
      ! detscale & detpiv - used for stats
      !

      pivsiz0: if (pivsiz.eq.0) then
         !print *, "pivsize0 dropout", ubest1, umin
         ! No pivot found in search of all available columns
         ! Since all the columns have been updated, revise m to q+1
         m = q+1
         kkm = kkq1
         rm = q
      end if pivsiz0

      pivsizes: if (pivsiz.eq.1) then
         !print *, "1x1 pivot on ", m
         nzero = 0 ! pivot found, reset recent zero count
         ! Swap columns if m not q+1
         if (q+1.ne.m) &
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm,q+1,m)
         ! We store D**-1.
         !print *, "1x1 pivot ", a(kkq1), "(col ", m, ")"
         d(2*q+1) = one/d11
         d(2*q+2) = zero
         if (a(kkq1).lt.0) stats%num_neg = stats%num_neg + 1
         ! Store L in A and LD in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = d11
         a(kkq1) = 1/d11
         !do i = 1, n-q-1
         !   buf(kq1+i) = a(kkq1+i)
         !   a(kkq1+i) = d(2*q+1)*a(kkq1+i)
         !end do
         pivval = a(kkq1) !d(2*q+1)
         jj = kb1+1
         do ii = kq1+1, kq1+p-q-1
            buf(jj) = a(ii); jj = jj + 1
            a(ii) = pivval*a(ii)
         end do
         do ii = kq1+p-q, kq1+n-q-1
            a(ii) = pivval*a(ii)
         end do
         ! Update columns q+2 to p
         kkj = kkq1
         do j = q+2, p
            kkj = kkj + lda + 1
            ii = kkq1+j-q-1
            alpha = buf(kb1+j-q-1)
            if(alpha.eq.zero) cycle
            !do k = 0, n-j
            !   a(kkj+k) = a(kkj+k) - alpha * a(ii+k)
            !end do
            do k = kkj, kkj+n-j-1, 2
               a(k)   = a(k)   - alpha * a(ii)
               a(k+1) = a(k+1) - alpha * a(ii+1)
               ii = ii + 2
            end do
            if(mod(n-j, 2).eq.0) then
               a(kkj+n-j) = a(kkj+n-j) - alpha * a(kkq1+j-q-1+n-j)
            end if
         end do
         ! Update q and kkq1
         kkq1 = kkq1 + lda + 1
         q = q+1

      else if (pivsiz.eq.2) then pivsizes
         !print *, "2x2 pivot on ", t, m
         nzero = 0 ! pivot found, reset recent zero count
         ! Swap columns unless t.eq.q+1 and m.eq.q+2
         if (q+2.ne.m) then
            if (q+1.ne.t) &
               call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm,q+1,t)
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm,q+2,m)
         end if
         ! We store D**-1
         k = kkq1 + lda + 1
         d(2*q+1) = (a(k)*detscale)/detpiv
         d(2*q+3) = (a(kkq1)*detscale)/detpiv
         d(2*q+2) = -(a(kkq1+1)*detscale)/detpiv
         d(2*q+4) = zero
         !print *, "2x2 pivot ", d(2*q+1:2*q+3)
         !print *, "   (cols ", t, m, ")"
         !print *, "   [vals ", a(k), a(kkq1), a(kkq1+1), &
         !   "dp", detscale, detpiv, "]"
         ! Update info
         if (detpiv<zero) then
            stats%num_neg = stats%num_neg + 1
         else if (a(kkq1)+a(k)<0) then
            stats%num_neg = stats%num_neg + 2
         end if
         stats%num_two = stats%num_two + 1
         ! Store L in A and LD or conjg(LD) in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = a(kkq1)
         buf(kb1+1) = a(kkq1+1)
         buf(kb1+ldbuf) = a(kkq1+1)
         buf(kb1+ldbuf+1) = a(k)
         a(kkq1) = 1.0_wp
         a(k) = 1.0_wp
         a(kkq1+1) = 0.0_wp
         do i = 2, p-q-1
            buf(kb1+i      ) = a(kkq1+i)
            buf(kb1+i+ldbuf) = a(k+i-1)
            a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
            a(k+i-1) = d(2*q+3)*a(k+i-1) + d(2*q+2)*buf(kb1+i)
         end do
         do i = p-q, n-q-1
            aval = a(kkq1+i)
            a(kkq1+i) = d(2*q+1)*a(kkq1+i) + d(2*q+2)*a(k+i-1)
            a(k+i-1) = d(2*q+3)*a(k+i-1) + d(2*q+2)*aval
         end do
         ! Update columns q+3 to p
         kkj = k + lda + 1
         j = max(q+3,1)
         if (p>=j) then
            call dgemm('n','t',n-j+1,p-j+1,2,-one,a(kkq1+j-q-1), &
               lda-((q+1)/nb)*nb,buf(kb1+j-q-1),ldbuf,one,a(kkj),lda)
         end if
         ! Update q and kkq1
         kkq1 = k + lda + 1
         q = q + 2

      else if (pivsiz.eq.-1) then pivsizes
         !print *, "ZERO pivot"
         nzero = nzero + 1 ! count zero pivot
         ! Handle a row that is zero
         if (q+1.ne.m) &
            call ma64_swap(n,nb,q,1,rm,a,lda,aleft,aln,buf,ldbuf,perm,q+1,m)
         !print *, "fail piv"
         d(2*q+1) = zero
         d(2*q+2) = zero
         if(.not.control%action) then
            ! Error on singular matrix
            stats%flag = MA97_ERROR_SINGULAR
            return
         end if
         stats%num_zero = stats%num_zero + 1
         ! Store L in A and -LD in buf
         kq1 = q+1+(q+0_long)*lda
         kb1 = q+1+(q+0_long)*ldbuf
         buf(kb1) = zero
         a(kkq1) = zero
         do i = 1, p-q-1
            buf(kb1+i) = zero
            a(kkq1+i) = zero
         end do
         do i = p-q, n-q-1
            a(kkq1+i) = zero
         end do
         kkq1 = kkq1 + lda + 1
         q = q+1
      end if pivsizes

      if (q .eq. p .or. pivsiz.eq.0) exit pivot
   end do pivot

   !print *, "exit d = ", d(1:2*q)
   !print *, "exit a = "
   !print "(5es12.4)", a
   !print *, "delayed ", p-q, "pivots"

end subroutine factor_solve_block_small

!*******************************
!
! performs factorization of trapezoidal matrix.
! based on hsl_ma64 version 6.0.0 (ported 18th January 2011)
!
subroutine factor_solve_block1(n, a, d, q, control, stats)
   integer, intent(in) :: n ! number of rows in trapeziodal matrix
   real(wp), intent(inout) :: a(*) ! holds trapezoidal matrix 
      ! to be factorized. Holds the factorization on exit.
      ! Delayed columns are permuted to the final columns.
   real(wp), intent(inout) :: d(2)
   integer, intent(out) :: q ! number of pivots chosen
   type(ma97_control), intent(in) :: control
   type(thread_stats), intent(inout) :: stats

   !     .. Local Scalars ..
   integer i ! Row index
   integer(long) k  ! Position in a
   real(wp) pivval ! temporary variable for storing pivotal value
   real(wp) u ! Relative pivot threshold
   real(wp) urel ! Relative pivot value

   u = min(max(control%u,zero),one)
   q = 0

   ! Check if diagonal is zero
   if (abs(a(1)).le.control%small) then
      ! Diagonal entry is effectively zero
      ! Check to see if entire column is zero

      ! Check entries in column. Stop as soon as we encounter a non-zero value
      do k = 2, n
         if (abs(a(k)).gt.control%small) then
            if (abs(a(k)) > huge(zero)) stats%flag = MA97_ERROR_INFINITY
            return ! Delay column
         end if
      end do

      ! Zero column
      ! All entries of the column are small
      if(.not.control%action) then
         ! Error on singular matrix
         stats%flag = MA97_ERROR_SINGULAR
         return
      end if
      stats%num_zero = stats%num_zero + 1
      ! Store L in A
      d(1:2) = zero
      do i = 1, n
         a(i) = zero
      end do
      q = q+1
      return
   end if

   ! Diagonal is non-zero
   ! Check entries in column. Stop as soon as any indicates bad pivot
   urel = abs(a(1)) / u
   do k = 2, n
      if (urel.le.abs(a(k))) return
   end do

   ! Good 1x1 pivot
   if (a(1).lt.0) stats%num_neg = stats%num_neg + 1
   ! Store L in A
   d(1) = 1/a(1)
   d(2) = zero
   a(1) = 1/a(1)
   pivval = a(1) !d(2*q+1)
   do i = 2, n
      a(i) = pivval*a(i)
   end do
   ! Update q
   q = q+1

end subroutine factor_solve_block1

!*******************************
!
! Swap columns. It may be assumed that 1<=j1<j2<=p.
! This was based on a routine within ma64.
!
subroutine ma64_swap(n,nb,q,q1,rm,a,lda,aleft,aln,buf,ldbuf,perm,j1,j2)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   integer, intent(in) :: q
   integer, intent(in) :: q1
   integer, intent(in) :: rm
   real(wp), intent(inout) :: a(*)
   integer, intent(in) :: lda
   real(wp), intent(inout) :: aleft(*)
   integer, intent(in) :: aln
   real(wp), intent(inout) :: buf(*)
   integer, intent(in) :: ldbuf
   integer, intent(inout) :: perm(*)
   integer, intent(in) :: j1,j2

   integer(long) d1,d2 ! Positions of the diagonal entries for columns j1,j2
   integer j ! Index of the leading column of the block
   integer jb ! Size of block
   integer(long) k1 ! Start of current part of column j1
   integer(long) k2 ! Start of current part of column j2
   integer l !
   real(wp) temp

   !print *, "swapping ", j1, j2

   j = perm(j1)
   perm(j1) = perm(j2)
   perm(j2) = j

   ! Swap columns of buffer
   call dswap(q-q1+1-rm,buf(j1+rm*ldbuf),ldbuf,buf(j2+rm*ldbuf),ldbuf)
   !call dswap(q-q1+1,buf(j1),n,buf(j2),n)

   k1 = j1
   k2 = j2

   ! Swap rows
   ! i.e. a(j1, 1:j1-1) <-> a(j2, 1:j1-1)
   l = j1-1
   if (l>0) call dswap(l,a(k1),lda,a(k2),lda)
   if (aln>0) call dswap(aln,aleft(j1),lda,aleft(j2),lda)

   ! calulate values for next phase
   d1 = k1 + l*int(lda,long)
   k1 = d1 + 1
   k2 = k2 + l*int(lda,long) + lda

   ! Swap columns with rows
   ! i.e. a(j1+1:j2-1, j1) <-> a(j2, j1+1:j2-1)
   jb = j2-j1-1
   l = min(jb,nb-j1,j2-1)
   if (l>0.and.j2.ge.1) call dswap(l,a(k1),1,a(k2),lda)

   ! set up values for next phase
   d2 = k2 + l*int(lda,long)
   k2 = d2 + 1
   k1 = k1 + l + 1

   ! Swap the diagonals
   ! i.e. a(j1,j1) <-> a(j2,j2)
   temp = a(d1)
   a(d1) = a(d2)
   a(d2) = temp

   ! Swap columns
   ! i.e. a(j2+1:n,j1) <-> a(j2+1:n,j2)
   if (n>j2) call dswap(n-j2,a(k1),1,a(k2),1)

end subroutine ma64_swap

!*******************************
!
! Update columns jl:jr for pivots pl:pr
! This was based on a routine within ma64.
!
subroutine ma64_update(n,nb,kkq,a,lda,buf,ldbuf,jl,jr,pl,pr)
   integer, intent(in) :: n
   integer, intent(in) :: nb
   integer(long), intent(in) :: kkq
   real(wp), intent(inout) :: a(*)
   integer, intent(in) :: lda
   real(wp), intent(in) :: buf(*)
   integer, intent(in) :: ldbuf
   integer, intent(in) :: jl,jr,pl,pr

   ! Local variables
   integer i
   integer j ! column index
   integer j2
   integer p
   integer(long) kkj ! Position of diagonal j

   integer, parameter :: bs=64 ! block size for update

   j = max(jl,1)
   j2 = min(nb,jr)

   if (n.lt.j .or. j2.lt.j .or. pl.gt.pr) return

   do i = 1, j2-j+1, bs
      kkj = 1 + (j-1+i-1)*(lda+1_long)
      p = min(bs, j2-j+1-i+1)
      call dgemm ('n', 't', n-j+1-i+1, p, pr-pl+1, -one, &
         a(kkq+j-1+lda*(pl-1)+i-1), lda, buf(j+ldbuf*(pl-1)+i-1), ldbuf, &
         one, a(kkj), lda)
   end do

end subroutine ma64_update

!*******************************
! Equivalent to call:
! call dgemv('NoTrans', m, n, -one, a, lda, b, ldb, one, c, 1)
! However for small updates exploit explicit loops as this is faster
! (calls dgemv itself for large updates)
subroutine mydgemv(m, n, a, lda, b, ldb, c)
   integer, intent(in) :: m
   integer, intent(in) :: n
   real(wp), dimension(*), intent(in) :: a
   integer, intent(in) :: lda
   real(wp), dimension(*), intent(in) :: b
   integer, intent(in) :: ldb
   real(wp), dimension(*), intent(inout) :: c

   integer :: i
   integer :: ipa, ipb
   real(wp) :: b1, b2, b3, b4, val

   select case(n)
   case(1)
      !call dgemv('NoTrans', m, 1, -one, a, lda, b, ldb, one, c, 1)
      b1 = b(1)
      c(1:m) = c(1:m) - b1 * a(1:m)
   case(2)
      !call dgemv('NoTrans', m, 2, -one, a, lda, b, ldb, one, c, 1)
      b1 = b(1)
      b2 = b(1+ldb)
      c(1:m) = c(1:m) - b1 * a(1:m) - b2 * a(lda+1:lda+m)
   case default
      select case(m)
      case(1)
         b1 = 0
         ipa = 1; ipb = 1
         do i = 1, n
            b1 = b1 + a(ipa) * b(ipb)
            ipa = ipa + lda
            ipb = ipb + ldb
         end do
         c(1) = c(1) - b1
      case(2)
         b1 = 0; b2 = 0
         ipa = 1; ipb = 1
         do i = 1, n
            val = b(ipb)
            b1 = b1 + a(ipa)   * val
            b2 = b2 + a(ipa+1) * val
            ipa = ipa + lda
            ipb = ipb + ldb
         end do
         c(1) = c(1) - b1
         c(2) = c(2) - b2
      case(3)
         b1 = 0; b2 = 0; b3 = 0
         ipa = 1; ipb = 1
         do i = 1, n
            val = b(ipb)
            b1 = b1 + a(ipa)   * val
            b2 = b2 + a(ipa+1) * val
            b3 = b3 + a(ipa+2) * val
            ipa = ipa + lda
            ipb = ipb + ldb
         end do
         c(1) = c(1) - b1
         c(2) = c(2) - b2
         c(3) = c(3) - b3
      case(4)
         b1 = 0; b2 = 0; b3 = 0; b4 = 0
         ipa = 1; ipb = 1
         do i = 1, n
            val = b(ipb)
            b1 = b1 + a(ipa)   * val
            b2 = b2 + a(ipa+1) * val
            b3 = b3 + a(ipa+2) * val
            b4 = b4 + a(ipa+3) * val
            ipa = ipa + lda
            ipb = ipb + ldb
         end do
         c(1) = c(1) - b1
         c(2) = c(2) - b2
         c(3) = c(3) - b3
         c(4) = c(4) - b4
      case default
         call dgemv('NoTrans', m, n, -one, a, lda, b, ldb, one, c, 1)
      end select
   end select
   
end subroutine mydgemv

!*************************************************
!
! Allocates memory on supplied stack
! ptr is set to point to a chunk of size len
! Note: no locks required as local to a single thread
!
subroutine stack_alloc(stack, ptr, len, st)
   type(stack_mem_type), pointer, intent(inout) :: stack
   real(wp), dimension(:), pointer, intent(out) :: ptr
   integer(long), intent(in) :: len
   integer, intent(out) :: st

   type(stack_mem_type), pointer :: stack_ptr

   st = 0

   ! If stack is null, allocate it
   if(.not.associated(stack)) then
      allocate(stack, stat=st)
      if(st.ne.0) return
      stack%mem_size = max(len, BLK_SZ)
      allocate(stack%mem(stack%mem_size), stat=st)
      if(st.ne.0) return
   end if

   ! Check if there is room on this page for the allocation
   if(stack%head+len.gt.stack%mem_size) then
      ! If insufficient space, allocate a new page at top of stack
      allocate(stack_ptr, stat=st)
      if(st.ne.0) return
      stack_ptr%below => stack
      stack => stack_ptr
      stack%mem_size = max(len, BLK_SZ)
      allocate(stack%mem(stack%mem_size), stat=st)
      if(st.ne.0) return
   end if

   ! Set pointer to next available block and increment head
   ptr => stack%mem(stack%head+1:stack%head+len)
   stack%head = stack%head + len
end subroutine stack_alloc

!*************************************************
!
! This releases the top len bytes of the stack. Frees should be used in
! reverse order to allocates (ie last on first off).
! If a stack page becomes empty, delete it.
! Note: no locks required as local to a single thread
!
subroutine stack_free(stack, len)
   type(stack_mem_type), pointer, intent(inout) :: stack
   integer(long), intent(in) :: len

   type(stack_mem_type), pointer :: stack_ptr

   stack%head = stack%head - len

   ! If a stack page is empty, delete it and update stack to next page down
   ! (which may be null)
   if(stack%head.eq.0) then
      stack_ptr => stack
      stack => stack%below
      deallocate(stack_ptr%mem)
      deallocate(stack_ptr)
   end if
end subroutine stack_free

!*************************************************
!
! Initialize memory storage of current unallocated page with specified sizes
! (minimum of BLK_SZ)
! MUST be called prior to smalloc() to set up alloc%lock
!
subroutine smalloc_setup(alloc, int_hint, real_hint, st)
   type(smalloc_type), intent(inout) :: alloc
   integer(long), intent(in) :: int_hint
   integer(long), intent(in) :: real_hint
   integer, intent(out) :: st

!$ call omp_init_lock(alloc%lock)

   st = 0
   if(.not.allocated(alloc%imem)) then
      alloc%imem_size = max(BLK_SZ+0_long,int_hint)
      allocate(alloc%imem(alloc%imem_size), stat=st)
   end if
   if(.not.allocated(alloc%rmem)) then
      alloc%rmem_size = max(BLK_SZ+0_long,real_hint)
      allocate(alloc%rmem(alloc%rmem_size), stat=st)
   end if
end subroutine smalloc_setup

!*************************************************
!
! Grab some real memory of the specified size and return a pointer to it
! Preference given to same page as last allocation.
! If this is not possible then check if we fit on any existing pages
! If not make a new page
! Note: as may be accessed by multiple threads, we have a lock
!
subroutine real_alloc(alloc_in, ptr, len, srcptr, srchead, st)
   type(smalloc_type), target, intent(inout) :: alloc_in
   real(wp), dimension(:), pointer, intent(out) :: ptr
   integer(long), intent(in) :: len
   type(smalloc_type), pointer, intent(out) :: srcptr
   integer(long), intent(out) :: srchead
   integer, intent(out) :: st

   type(smalloc_type), pointer :: alloc

   st = 0
   if(len.lt.0) return

!$ call omp_set_lock(alloc_in%lock)

   ! First try same page as last alloc
   if(associated(alloc_in%top_real)) then
      alloc => alloc_in%top_real
      if(alloc%rhead+len.le.alloc%rmem_size) then
         ! Sufficient space, allocate in this block
         ptr => alloc%rmem(alloc%rhead+1:alloc%rhead+len)
         srcptr => alloc
         srchead = alloc%rhead+1
         alloc%rhead = alloc%rhead + len
!$       call omp_unset_lock(alloc_in%lock)
         return
      end if
   end if

   ! Else check all pages, if we reach the end create a new one of sufficient
   ! size to allocate the pointer
   alloc => alloc_in
   do
      if(.not.allocated(alloc%rmem)) then
         alloc%rmem_size = max(len,BLK_SZ)
         allocate(alloc%rmem(alloc%rmem_size), stat=st)
         if(st.ne.0) then
!$          call omp_unset_lock(alloc_in%lock)
            return
         end if
      end if
      if(alloc%rhead+len.le.alloc%rmem_size) then
         ! Sufficient space, allocate in this block
         ptr => alloc%rmem(alloc%rhead+1:alloc%rhead+len)
         srcptr => alloc
         srchead = alloc%rhead+1
         alloc%rhead = alloc%rhead + len
         alloc_in%top_real => alloc
         exit
      end if
      ! Insufficent space, move to next block
      if(.not.associated(alloc%next_alloc)) then
         allocate(alloc%next_alloc, stat=st)
         if(st.ne.0) then
!$          call omp_unset_lock(alloc_in%lock)
            return
         end if
      end if
      alloc => alloc%next_alloc
   end do

!$ call omp_unset_lock(alloc_in%lock)
end subroutine real_alloc

!*************************************************
!
! Grab some integer memory of the specified size and return a pointer to it
! Preference given to same page as last allocation.
! If this is not possible then check if we fit on any existing pages
! If not make a new page
! Note: as may be accessed by multiple threads, we have a lock
!
subroutine int_alloc(alloc_in, ptr, len, srcptr, srchead, st)
   type(smalloc_type), target, intent(inout) :: alloc_in
   integer, dimension(:), pointer, intent(out) :: ptr
   integer(long), intent(in) :: len
   type(smalloc_type), pointer, intent(out) :: srcptr
   integer(long), intent(out) :: srchead
   integer, intent(out) :: st

   type(smalloc_type), pointer :: alloc

   st = 0
   if(len.lt.0) return

!$ call omp_set_lock(alloc_in%lock)

   ! First try same page as last alloc
   if(associated(alloc_in%top_int)) then
      alloc => alloc_in%top_int
      if(alloc%ihead+len.le.alloc%imem_size) then
         ! Sufficient space, allocate in this block
         ptr => alloc%imem(alloc%ihead+1:alloc%ihead+len)
         srcptr => alloc
         srchead = alloc%ihead+1
         alloc%ihead = alloc%ihead + len
!$       call omp_unset_lock(alloc_in%lock)
         return
      end if
   end if

   ! Else check all pages, if we reach the end create a new one of sufficient
   ! size to allocate the pointer
   alloc => alloc_in
   do
      if(.not.allocated(alloc%imem)) then
         alloc%imem_size = max(len,BLK_SZ)
         allocate(alloc%imem(alloc%imem_size),stat=st)
         if(st.ne.0) then
!$          call omp_unset_lock(alloc_in%lock)
            return
         end if
      end if
      if(alloc%ihead+len.le.size(alloc%imem)) then
         ! Sufficient space, allocate in this block
         ptr => alloc%imem(alloc%ihead+1:alloc%ihead+len)
         srcptr => alloc
         srchead = alloc%ihead+1
         alloc%ihead = alloc%ihead + len
         alloc_in%top_int => alloc
         exit
      end if
      ! Insufficent space, move to next block
      if(.not.associated(alloc%next_alloc)) then
         allocate(alloc%next_alloc, stat=st)
         if(st.ne.0) then
!$          call omp_unset_lock(alloc_in%lock)
            return
         end if
      end if
      alloc => alloc%next_alloc
   end do

!$ call omp_unset_lock(alloc_in%lock)
end subroutine int_alloc

!*************************************************
!
! Free all memory associated with the specified alloc linked list
! (both real and integer)
!
subroutine smfreeall(alloc_in)
   type(smalloc_type), intent(inout) :: alloc_in

   type(smalloc_type), pointer :: alloc, alloc2
   integer :: st

   !integer :: cnt
   !integer :: wasted

!$ call omp_destroy_lock(alloc_in%lock)

   deallocate(alloc_in%rmem, stat=st)
   alloc_in%rhead = 0
   deallocate(alloc_in%imem, stat=st)
   alloc_in%ihead = 0
   alloc => alloc_in%next_alloc
   nullify(alloc_in%next_alloc, alloc_in%top_real, alloc_in%top_int)
   !cnt = 1
   !wasted = wasted + BLK_SZ - alloc_in%rhead

   do while(associated(alloc))
      !cnt = cnt + 1
      !wasted = wasted + BLK_SZ - alloc%rhead
      deallocate(alloc%rmem, stat=st)
      deallocate(alloc%imem, stat=st)
      alloc2 => alloc%next_alloc
      deallocate(alloc)
      alloc => alloc2
   end do
   !print *, "Freed", cnt, "pages", wasted/real(BLK_SZ*cnt), "wasted."
end subroutine smfreeall

!*************************************************
!
! routine to print errors and warnings
!
subroutine ma97_print_flag(context,nout,iflag,st)
   integer, intent(in) :: iflag, nout
   integer, intent(in), optional :: st
   character (len=*), optional, intent(in) :: context

   if (nout < 0) return
   if (iflag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(context),&
         '. Error flag = ', iflag
   else
      write (nout,'(/3a,i3)') ' Warning from ',trim(context),&
         '. Warning flag = ', iflag
   end if

   ! Errors
   select case(iflag)
   case(MA97_ERROR_CALL_SEQUENCE)
      write (nout,'(a)') ' Error in sequence of calls.'

   case(MA97_ERROR_A_N_OOR)
      write (nout,'(a)') ' n or ne is out of range (or has changed)'

   case(MA97_ERROR_A_PTR)
      write (nout,'(a)') ' Error in ptr'

   case(MA97_ERROR_A_ALL_OOR)
      write (nout,'(a)') ' All entries in a column out-of-range (ma97_analyse)'
      write (nout,'(a)') ' or all entries out-of-range (ma97_analyse_coord)'

   case(MA97_ERROR_MATRIX_TYPE)
       write (nout,'(a)') ' matrix_type is out of range or has changed'

   ! Below error is only present in complex version
   !case(MA97_ERROR_IMAG_DIAGONAL)
   !   write (nout,'(a)') ' one or more diagonal entries is not real'

   case(MA97_ERROR_SINGULAR)
      write (nout,'(a)') ' Matrix found to be singular'

   case(MA97_ERROR_NOT_POS_DEF)
       write (nout,'(a)') ' Matrix is not positive-definite'

   case(MA97_ERROR_INFINITY)
      write (nout,'(a)') ' IEEE infinities detected. Factorization terminated.'

   case(MA97_ERROR_PTR_ROW)
      write (nout,'(a)') ' ptr and row should be present'

   case(MA97_ERROR_ORDER)
      write (nout,'(a/a)') &
         ' Either control%ordering out of range or error in user-supplied  &
         &elimination order'

   case(MA97_ERROR_X_SIZE)
      write (nout,'(a)') ' Error in size of x or nrhs'

   case(MA97_ERROR_JOB_OOR)
      write (nout,'(a,i10)') ' job out of range.'

   case(MA97_ERROR_NOT_LLT)
      write (nout,'(a)')&
         ' Not a LL^T factorization of a positive-definite matrix'

   case(MA97_ERROR_NOT_LDLT)
      write (nout,'(a)')&
         ' Not a LDL^T factorization of an indefinite matrix'

   case(MA97_ERROR_ALLOCATION)
      if (present(st)) then
         write (nout,'(a,i6)') ' Allocation error. stat parameter = ', st
      else
         write (nout,'(a)') ' Allocation error'
      end if

   case(MA97_ERROR_NO_METIS)
      write (nout,'(a)') &
         ' METIS requested but not available'

   case(MA97_ERROR_MC68)
      write (nout,'(a)') &
         ' Unexpected error return from HSL_MC68 (called by ma97_analyse)'

   case(MA97_ERROR_MC77)
      write (nout,'(a)') &
         ' Unexpected error return from MC77 (called by ma97_factor)'

   case(MA97_ERROR_VAL)
      write (nout,'(a)') &
         ' Optional argument val not present when expected'

   case(MA97_ERROR_NO_SAVED_SCALING)
      write (nout,'(a)') &
         ' Requested use of scaling from matching-based &
         &ordering but matching-based ordering not used.'

   case(MA97_ERROR_NBI_OOR)
      write (nout,'(a)') ' nbi out of range.'

   ! Warnings
   case(MA97_WARNING_IDX_OOR)
      write (nout,'(a)') ' out-of-range indices detected'

   case(MA97_WARNING_DUP_IDX)
      write (nout,'(a)') ' duplicate entries detected'

   case(MA97_WARNING_DUP_AND_OOR)
      write (nout,'(a)') &
         ' out-of-range indices detected and duplicate entries detected'

   case(MA97_WARNING_MISSING_DIAGONAL)
      write (nout,'(a)') ' one or more diagonal entries is missing'

   case(MA97_WARNING_MISS_DIAG_OORDUP)
      write (nout,'(a)') ' one or more diagonal entries is missing and'
      write (nout,'(a)') ' out-of-range and/or duplicate entries detected'

   case(MA97_WARNING_ANAL_SINGULAR)
      write (nout,'(a)') ' Matrix found to be structually singular'

   case(MA97_WARNING_FACT_SINGULAR)
      write (nout,'(a)') ' Matrix found to be singular'

   case(MA97_WARNING_MATCH_ORD_NO_SCALE)
      write (nout,'(a)') &
         ' Matching-based ordering used but associated scaling ignored'

   case default
      write (nout,'(a)') ' HSL_MA97 Internal Error '

   end select

end subroutine ma97_print_flag

!*************************************************
!
! This subroutine converts a node number X to a string "nX" without space
! between the n and the X.
! It is used to output the graph in a nice fashion if dograph=.true.
!
!%%subroutine node_name(str, num)
!%%   character(len=8), intent(out) :: str
!%%   integer, intent(in) :: num
!%%
!%%   str(1:8) = "        "
!%%   select case(num)
!%%   case(0:9)
!%%      write(str, "(a1,i1)") "n", num
!%%   case(10:99)
!%%      write(str, "(a1,i2)") "n", num
!%%   case(100:999)
!%%      write(str, "(a1,i3)") "n", num
!%%   case(1000:9999)
!%%      write(str, "(a1,i4)") "n", num
!%%   case(10000:99999)
!%%      write(str, "(a1,i5)") "n", num
!%%   case default
!%%      print *, "node name too large"
!%%      stop
!%%   end select
!%%end subroutine node_name

!
! This subroutine can be used to compare to fkeeps
! Used for debugging between different parallel runs to check answer is the same
!
!subroutine cmp_fkeep(akeep, fkeep1, fkeep2)
!   type(ma97_akeep), intent(in) :: akeep
!   type(ma97_fkeep), intent(in) :: fkeep1
!   type(ma97_fkeep), intent(in) :: fkeep2
!
!   integer :: nd
!   integer :: node
!   integer :: blkm
!   integer :: blkn
!   integer :: nelim
!   logical :: diff
!   integer(long) :: ip1, ip2
!   integer :: i 
!
!   real(wp), dimension(:), pointer :: lcol1, lcol2
!
!   print *, "Checking nodes:"
!   do node = 1, akeep%nnodes
!      diff = .false.
!      if(fkeep1%nodes(node)%nelim .ne. fkeep2%nodes(node)%nelim) then
!         print *, "Node ", node, "has differing nelim:", &
!            fkeep1%nodes(node)%nelim, fkeep2%nodes(node)%nelim
!         diff = .true.
!      end if
!      if(fkeep1%nodes(node)%ndelay .ne. fkeep2%nodes(node)%ndelay) then
!         print *, "Node ", node, "has differing ndelay:", &
!            fkeep1%nodes(node)%ndelay, fkeep2%nodes(node)%ndelay
!         diff = .true.
!      end if
!      if(any(fkeep1%nodes(node)%perm .ne. fkeep2%nodes(node)%perm)) then
!         print *, "Node ", node, "has differing perm:"
!         print *, "Perm 1:", fkeep1%nodes(node)%perm
!         print *, "Perm 2:", fkeep2%nodes(node)%perm
!         diff = .true.
!      end if
!      if(diff) then
!         print *, "   (not checking lcol)"
!      else
!         nd = fkeep1%nodes(node)%ndelay
!         nelim = fkeep1%nodes(node)%nelim
!         blkn = akeep%sptr(node+1) - akeep%sptr(node) + nd
!         blkm = int(akeep%rptr(node+1) - akeep%rptr(node)) + nd
!         ! lcol
!         if(any(fkeep1%nodes(node)%lcol .ne. fkeep2%nodes(node)%lcol)) then
!            print *, "========"
!            print *, "Node ", node, "has differing lcol."
!            print *, "Is ", blkm, "x", blkn, "with ", nd, "delays in"
!            print *, "nelim = ", nelim
!            print *, "parent = ", akeep%sparent(node)
!            lcol1 => fkeep1%nodes(node)%lcol
!            lcol2 => fkeep2%nodes(node)%lcol
!            do i = 1, nelim
!               ip1 = (i-1)*blkm + i
!               ip2 = i*blkm
!               if(any(lcol1(ip1:ip2).ne.lcol2(ip1:ip2))) then
!                  print *, "   col ", i, "differs:"
!                  print *, "      1:", lcol1(ip1:ip2)
!                  print *, "      2:", lcol2(ip1:ip2)
!                  print *, "   diff:", abs(lcol1(ip1:ip2)-lcol2(ip1:ip2))
!                  print *, "   [last entry ip = ", ip2, "]"
!               end if
!            end do
!         end if
!      end if
!   end do
!   print *, "== End of node diff =="
!end subroutine cmp_fkeep

pure integer function ma97_get_n_double(akeep)
   type(ma97_akeep), intent(in) :: akeep

   ma97_get_n_double = akeep%n
end function ma97_get_n_double

pure integer function ma97_get_nz_double(akeep)
   type(ma97_akeep), intent(in) :: akeep

   ma97_get_nz_double = akeep%ne
end function ma97_get_nz_double

end module hsl_ma97_double
