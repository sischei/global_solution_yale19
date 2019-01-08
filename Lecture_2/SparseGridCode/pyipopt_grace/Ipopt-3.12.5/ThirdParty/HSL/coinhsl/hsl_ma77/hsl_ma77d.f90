! COPYRIGHT (c) 2006 Council for the Central Laboratory
!               of the Research Councils
! Original date 5 September 2006. Version 1.0.0.

! IMPORTANT : Version 1.0.0 for HSL 2007 was ONLY for the positive definite
!             Version 2.0.0 is for positive definite and indefinite matrices
!             Version 3.0.0 allows for running on 64-bit architectures
!             Version 4.0.0 option to scale added
!             Version 5.0.0 static pivoting modified (following 
!                           change to hsl_ma64)

! Written by:  John Reid and Jennifer Scott

! Version 5.8.0
! For full version history see ChangeLog

 ! Some general remarks:
!  Elements are linked by earliest variable in pivot sequence.
!  Supervariables used so that when a variable is eliminated in
!  the row-entry code, variables with identical patterns are too.
!  The user decides on entry to ma77_factor whether the problem is to be
!  treated as pos.def. or indefinite.
!  We allow 2x2 pivots to be input in the indefinite version.
!  Frontal matrix is held in main memory.
!

!!!!!
! To generate the single precision version from the double version:
! (1) globally change _double to _single
! (2) change 0.0d0 to 0.0
! (3) globally change dcopy, daxpy  to scopy, saxpy, respectively.
! (4) globally change huge(0_short)/8 to huge(0_short)/4
!                 and huge(0_long)/8  to huge(0_long)/4.
!!!!!

module hsl_MA77_double

  use hsl_ma54_double
  use hsl_ma64_double
  use hsl_of01_double, of01_rdata => of01_data
  use hsl_of01_integer, of01_idata => of01_data

  implicit none

!*************************************************
! Parameters (all private)
! All integers are declared to be either short or long
  integer, parameter, private  :: wp = kind(0.0d0)
  integer, parameter, private  :: long = selected_int_kind(18)
  integer, parameter, private  :: short = kind(0)
  real (wp), parameter, private :: one = 1.0_wp
  real (wp), parameter, private :: zero = 0.0_wp
  integer(short), parameter, private :: nemin_default = 8
  integer(short), parameter, private :: lup1 = huge(0_short)/8 !(do NOT 
!    change this for 64-bit version)
  integer(short), parameter, private :: nb54_default = 150
  integer(short), parameter, private :: nb64_default = 120
  integer(short), parameter, private :: nbi_default = 40
!*************************************************

  interface MA77_open
      module procedure MA77_open_double
  end interface

  interface MA77_input_vars
      module procedure MA77_input_vars_double
  end interface

  interface MA77_input_reals
      module procedure MA77_input_reals_double
  end interface

  interface MA77_analyse
      module procedure MA77_analyse_double
  end interface

  interface MA77_factor
      module procedure MA77_factor_double
  end interface

  interface MA77_factor_solve
      module procedure MA77_factor_solve_double
  end interface

  interface MA77_solve
      module procedure MA77_solve_double
  end interface

  interface MA77_resid
      module procedure MA77_resid_double
  end interface

  interface MA77_scale
      module procedure MA77_scale_double
  end interface

  interface MA77_enquire_posdef
    module procedure MA77_enquire_posdef_double
  end interface

  interface MA77_enquire_indef
    module procedure MA77_enquire_indef_double
  end interface

  interface MA77_alter
    module procedure MA77_alter_double
  end interface

  interface MA77_restart
      module procedure MA77_restart_double
  end interface

  interface MA77_finalise
      module procedure MA77_finalise_double
  end interface

!****************************************************************************

  type MA77_control ! The scalar control of this type controls the action

    logical :: action = .true. ! pos_def = .false. only.
!         If true and the matrix is found to be
!         singular, computation continues with a warning.
!         Otherwise, terminates with error -11.

    integer(short) :: bits = 32 ! set to 64 for 64-bit architecture

    integer(short) :: buffer_lpage(2) = 2**12 ! Number of scalars held in each
!           page of the integer and real in-core buffers.

    integer(short) :: buffer_npage(2) = 1600 ! Number of pages in the integer
!           and real in-core buffers.

    integer(long) :: file_size = 2**21 ! Target filesize for each of01 file.

    integer(short) :: infnorm = 0 ! controls norm used by scaling. Uses
!           one norm if infnorm /= 0 and infinity norm if infnorm = 0.

    integer(short) :: maxit = 1 ! max number of iterations performed
!          by scaling algorithm

    integer(long) :: maxstore = 0_long ! Max. amount of storage (in Fortran
!           storage units) for in-core arrays that replace superfiles.

    real(wp) :: multiplier = 1.1 ! In the indefinite case, the predicted
!           frontsize is increased by max(1,multipler) at start of
!           factorize and whenever front found to be too small.

    integer(short) :: nb54 = 150 ! Block size for the hsl_ma54 kernel code.

    integer(short) :: nb64 = 120 ! Block size for the hsl_ma64 kernel code.

    integer(short) :: nbi = 40 ! Inner block size for the hsl_ma64 kernel code

    integer(short) :: nemin = nemin_default !
!           Min. number of eliminations at a tree
!           node for amalgamation not to be considered.

   integer(short) :: p = 4 ! Controls choice of p:
!           p=1: Split point at 1
!           p=2: Split point at 2
!           p=3: Split point at nc
!           p=4: Split point chosen by G & L'E

    integer(short) :: print_level = 0 ! Controls diagnostic printing.
!           Possible values are:
!                       < 0: no printing.
!                         0: error and warning messages only.
!                         1: as 0 plus basic diagnostic printing.
!                       > 1: as 1 plus some more detailed diagnostic messages.

    real(wp) :: small = 1e-20_wp ! Minimum pivot size (absolute value of a
!           pivot must be of size at least small to be accepted).

    real (wp) :: static = zero ! Controls static pivoting. Only used
!           in indefinite case. If static > 0, pivots that do not
!           satisfy the threshold criteria may be selected and small
!           pivots may be replaced by static. In this case, the
!           factorization may be inaccurate.

    integer(long) :: storage(3) = 0_long ! If positive, numbers of scalars
!           initially held in the arrays that replace superfiles.

    integer(long) :: storage_indef = 0 ! If positive, numbers of scalars
!           initially held in the arrays that replace superfile
!           that is used as workspace in indefinite case.

    real(wp) :: thresh = 0.5 ! Tolerance for stopping the iterations
!           of the iterative scaling algorithm

    integer(short) :: unit_diagnostics = 6 ! unit number for
!           diagnostic printing.
!           Printing is suppressed if unit_diagnostics  <  0.

    integer(short) :: unit_error = 6 ! unit number for error messages
!           Printing is suppressed if unit_error  <  0.

    integer(short) :: unit_warning = 6 ! unit number for warning messages
!           Printing is suppressed if unit_warning  <  0.

    real (wp) :: u = 0.01 ! Initial relative pivot threshold. Only used
!           in indefinite case.
!           This is a control parameter for HSL_MA64.

    real (wp) :: umin = 1.00 ! Minimum relative pivot threshold. Only used
!           in indefinite case. If umin < u, the pivot test may
!           be relaxed to prevent pivots from being delayed
!           (that is u may be replaced by a smaller value, which is
!           at least umin).
!           This is a control parameter for HSL_MA64.

  end type MA77_control

!****************************************************************************

  type MA77_info ! The scalar info of this type returns information to user.
    real (wp) :: detlog = 0.0 ! log. of absolute value of the determinant
    integer(short) :: detsign = 1  ! Sign of the determinant of the matrix.
    integer(short) :: flag = 0     ! Value zero after successful entry.
!            Possible error returns:
!              -1 Allocation error
!              -3 Error in the sequence of calls
!              -4 n < 0
!              -5 Error from of01 (error in Fortran inquire)
!              -6 Error from of01 (error in Fortran read)
!              -7 Error from of01 (error in Fortran open)
!              -8 Deallocation error
!              -9 MA77_input_vars has already been called for this element/row
!             -10 MA77_input_vars was not called for this element/row
!             -11 problem found to be singular or unexpectedly not pos. def
!             -12 Error from of01 (a filename with same name already exists)
!             -13 len(filename) (or len(restart_file)) is too large
!             -14 data for previous element/row incomplete
!             -15 Error from of01 (error in Fortran write)
!             -16 len(path) is too large or
!                 open unsuccessful for all elements of path.
!             -17 out-of-range and/or duplicated indices entered and
!                 user has not entered all the reals for the current elt/row
!                 in a single call to inelrs.
!             -18 nelt < 0
!             -19 Length of array reals is < 0
!             -20 job is out of range
!             -21 Error in user-supplied pivot order held in order.
!             -22 For one or more elements/rows, MA77_input_vars was called
!                 but no corresponding call was made to MA77_input_reals.
!             -23 control%buffer_lpage(:) < 1 or > control%file_size
!             -24 Error in size of x.
!             -25 Error in size of resid.
!             -26 control%buffer_npage(:) < 1
!             -27 control%maxstore out of range
!             -28 Files with names restart_file and/or filename(1:2)
!                 do not exist.
!             -29 pos_def = .true. but 2x2 pivots supplied to MA77_analyse
!             -30 Error return MA77_factor if front too large for allocation
!             -31 All the variables in an element/row are out-of-range
!             -32 too many reals input for current element/row
!             -33 number of entries in variable list < 0
!             -34 Returned by MA77_enquire_indef or MA77_alter
!                 if the call does not follow a successful call to
!                 MA77_factor with pos_def = .false.
!             -35 Returned by MA77_enquire_posdef
!                 if the call does not follow a successful call to
!                 MA77_factor with pos_def = .true.
!             -36 Returned by MA77_enquire_posdef,
!                 MA77_enquire_indef, or MA77_alter
!                 if there is an error in the size of array piv_order.
!             -37 Returned by MA77_enquire_posdef,
!                 MA77_enquire_indef, or MA77_alter
!                 if there is an error in the size of array d
!             -38 static < small and static /= 0.0
!             -39 length array scale too small, or scale present when
!                 should not be or absent when expected on call to ma77_solve
!             -40 IEEE infinities found in the reduced matrix
!             -41 Error in  Fortran close statement

!           Possible warnings:
!             +1 Out-of-range variable indices
!             +2 Duplicated variable indices
!             +4 Problem singular (control%action = .true.)
    integer(short) :: iostat = 0 ! iostat value on error
!          return -5, -6, -7, -15.
    integer(short) :: matrix_dup = 0 ! Number of duplicated entries.
    integer(short) :: matrix_rank = 0 ! Rank of factorized matrix.
    integer(short) :: matrix_outrange = 0 ! Number of out-of-range entries.
    integer(short) :: maxdepth = 0 ! Maximum depth of the tree.
    integer(short) :: maxfront = 0 ! Maximum front size.
    integer(long)  :: minstore = 0_long ! Amount of storage used in superfiles
    integer(short) :: ndelay = 0 ! Number of delayed eliminations.
    integer(long)  :: nfactor = 0_long ! Number of entries in the factor.
    integer(long)  :: nflops = 0_long ! Number of flops needed to calculate L.
    integer(short) :: niter = 0 ! Number of iterations of scaling algorithm.
    integer(short) :: nsup = 0 ! Number of supervariables. (Copy of keep%nsup).
    integer(short) :: num_neg = 0 ! Number of negative eigenvalues.
    integer(short) :: num_nothresh = 0 ! Number of pivots which did not
!           satisfy the threshold criteria.
    integer(short) :: num_perturbed = 0 ! holds number of pivots that were
!          replaced by control%static
    integer(short) :: ntwo = 0 ! Number of 2x2 pivots.
    integer(short) :: stat = 0 ! STAT value on error return -1 and -8.
    integer(short) :: index(1:4) = -1
!                  Holds info. on index of of01 superfiles used
!                  index(1) index of main integer superfile
!                  index(2) index of main real superfile
!                  index(3) index of real work superfile
!                  index(4) index of real superfile for delayed pivot cols
!                  index(i) < 0 if no data written to superfile.

    integer(long) :: nio_read(1:2) = 0_long ! On exit from a call to 
!          MA77_analyse, MA77_factor, and MA7_solve, holds the no. of integer 
!          and real records read by HSL_OF01 during the subroutine call
    integer(long) :: nio_write(1:2) = 0_long ! On exit from a call to 
!          MA77_analyse, MA77_factor, and MA7_solve, holds the no. of integer
!          and real records written by HSL_OF01 during the subroutine call
    integer(long) :: nwd_read(1:2) = 0_long ! On exit from a call to 
!          MA77_analyse, MA77_factor, and MA7_solve, holds the no. of integer
!          and real scalars read by HSL_OF01 during the subroutine call
    integer(long) :: nwd_write(1:2) = 0_long ! On exit from a call to 
!          MA77_analyse, MA77_factor, and MA7_solve, holds the no. of integer
!          and real scalars written by HSL_OF01 during the subroutine call

    integer(short) :: num_file(1:4) = 0 ! On exit from a call to MA77_finalise,
!          holds the number of secondary files that have been used
!          by the superfiles

    integer(long) :: storage(1:4) = 0_long ! Holds info. on storage used for
!          arrays that replace superfiles.
!          storage(1) ! no. of integers stored in main integer superfile
!                       (set to keep%ifree-1 at end of factorize)
!          storage(2) ! no. of reals stored in main real superfile
!                       (set to keep%rfree-1 at end of factorize)
!          storage(3) ! max. size of main real stack (keep%rtopmx)
!          storage(4) ! max. size of real stack for delayed pivot
!                       cols (keep%rtopdmx)
    integer(short) :: tree_nodes = 0 ! number of non-leaf nodes in tree
!           (including any that are discarded but not reused).
!           (=keep%tnode-keep%nelt)
    integer(short) :: unit_restart = -1 ! Holds unit number of restart file
    integer(short) :: unused = 0 ! Holds number of unused variables.
    real(wp) :: u = zero ! if num_perturbed = 0, on exit from factor
!          u holds threshold parameter that was used.

!     integer(short) :: wasted = 0
  end type MA77_info

!****************************************************************************

  type MA77_node ! Each element of the array tree of this type holds
!                data at a non-leaf node of the tree.
!    private ! 22.07.08 no longer private because hsl_ma79 needs to
!              have access to some of its components.
    integer(short), allocatable :: child(:) ! Names of the node's children.
! Note: children are all elements (either original elements or
!       generated elements: original rows are NOT included in the
!       tree structure).
    integer(short) :: nelim ! Number of variables eliminated at the node if
!          there are no delayed pivots.
  end type MA77_node

!****************************************************************************

! Structure of the main integer superfile is as follows:
! (a) For each row/element i input by the user in the
!     calls to MA77_input_vars we store
!     a list of length nvar = keep%size(i) of user-supplied variables.
!     If there are duplicates and/or out-of-range entries, we flag this
!     by setting keep%size(i) = -nvar, where nvar is the number of variables
!     in list AFTER duplicates/out-of-range entries squeezed out.
!     We write this list out to file, then the no. of entries
!     in the original user-supplied list and then a mapping
!     from the original list to the compressed list.
!     map(j) = 0 indicates jth variable out-of-range and to be ignored.
!     if map(j) = map(k) then  jth and kth variables are duplicates.
! (b) For each non-leaf node we store the list of variables
!     in elimination order,
!     with those that are eliminated at the node given first.
!     keep%ifile(i) holds the start of the data for node i.
! (c) In the indefinite case, during MA77_factor
!     we hold the actual lists
!     of the variables at each node (which will differ from
!     those stored during analyse if pivots are delayed).
!     keep%posint holds the position beyond the
!     end of the variable list for the last node in the tree
!     generated by MA77_analyse; this is the start of the lists
!     of integers for the nodes of the tree used during
!     the factorization. If more than one problem is
!     factorized, the code will overwrite from keep%posint.
!     Note: if the user has input 2x2 pivots, then negative signs
!     are used to flag this (amking it necessary to take absolute
!     values at non-leaf nodes).

! Structure of the main real superfile is as follows:
! (a) For each row/element i input by the user in the
!     calls to MA77_input_reals we store the real data supplied
!     by the user (nvar reals in row case and
!     (nvar*nvar+nvar)/2 in element case, where
!     nvar = keep%size(i)). keep%rfile(i) holds the position after
!     the last entry in the list for elt/equ i (it points
!     to the end of the list to allow user to input data using
!     more than one call)
! (b) keep%posfac holds the position beyond the
!     end of the reals for the last element/row
!     input by the user; this is the start of the reals on diagonal (block)
!     of factor followed by for the factors. In positive definite
!     case, we store n diagonal entries (square roots
!     of the pivots), in the indefinite case we leave
!     room for 2n entries (to allow for 2x2 pivots) before we write the
!     factor entries (actually store the inverse of D in indefinite case)
!     If more than one problem is
!     factorized, the code will overwrite from keep%posfac.

!****************************************************************************

  type MA77_keep ! The scalar keep of this type holds matrix data that is
!                  passed between subroutines
!    private ! 22.07.08 keep is no longer private because hsl_ma79 needs to
!              have access to some of its components.
    integer(long) :: dfree ! Last location in the main real superfile
!                            used for storage of diagonal entries.
    logical :: element_input ! .true. iff input is by elements.
    integer(long) :: file_size ! file_size for of01
    integer(short) :: flag = 0 ! holds warning flag (needed so that on reverse
!           communications calls, info%flag can be reset)
    integer(long) :: ifree ! First free location in the main integer superfile
    integer(short) :: index(1:4) = -1 ! Holds index of of01 superfiles used
!              index(1) index of main integer superfile
!              index(2) index of main real superfile
!              index(3) index of real work superfile
!              index(4) index of real superfile for delayed pivot columns
!              index(i) < 0 if no data written to file.
!   In pos. def. case, files are opened on index(4:5) but not actually used.
    integer(short) :: inelrs ! Number of reals required to complete the input
!           of the current element/row (used by MA77_input_reals).
    integer(short) :: inelrn ! index of the element/row that was last input.
    integer(short) :: lpage(2) ! lpage for of01
    integer(short) :: ltree ! Size of the array called tree.
    integer(long)  :: lup ! huge(0_short)/8 (32-bit) or huge(0_long)/8 (64-bit)
    integer(short) :: l1,l2 ! we use tree(keep%l1:keep%l2)
!               l1 is set to keep%nelt+1
!               l2 is set to keep%nelt+keep%ltree
!               If finalise is called with restart_file present,
!               l2 is set to keep%tnode (we allocate tree(keep%l1:keep%l2)
!               for restarting as the actual number of nodes in
!               tree is known by this stage)
    integer(short) :: matrix_dup ! Number of duplicated entries.
    integer(short) :: matrix_outrange ! Number of out-of-range entries.
    integer(short) :: maxelim ! Holds largest number of eliminations at a node
    integer(short) :: maxelim_actual ! Set in factorize to hold
!           the largest number of eliminations
!           actually performed at a node (indefinite case only)
    integer(long)  :: maxfa ! Holds the largest factor contribution that
!           is written to the main real superfile during call to MA77_factor.
!           We need this to allocate sufficient space for
!           reading in factor blocks during solve.
    integer(short) :: maxfront ! Holds the largest front size.
!           Set during call to MA77_analyse. Not altered during factor.
    integer(short) :: maxdepth ! Holds max tree depth.
    integer(short) :: maxfrontb ! Holds the largest front size
!           during factorize (=maxfront is pos def case)
    integer(short) :: maxlen   ! Holds the number of variables in
!           the largest node in the tree.
!           Set during call to MA77_analyse. Not altered during factor.
    integer(long) :: maxstore ! Max. amount of storage available
!           (in Fortran storage units) for in-core working
    integer(short) :: mvar ! Max. number of variables in an elt/row.
    integer(short) :: mmvar ! Max. number of variables in an elt/row after
!           duplicates and out-of-range entries removed.
    integer(long)  :: mx_ifree ! mx_ifree is
!           set to no. of integers written to integer superfile
!           (needed so that info%storage(1) is calculated correctly)
    integer(short) :: n ! Order of the system.
!           Set by user on the call to MA77_open.
    character(50)  :: name ! Procedure name (used when printing messages).
    integer(short) :: nb ! Copy of control%nb54 or control%nb64
!           (set on call to MA77_factor).
    integer(short) :: nbi ! Copy of control%nbi
!           (set on call to MA77_factor).
    integer(short) :: nelt ! Bound on the integers used to index elements in
!           the element-entry case or variables in the row-entry
!           case. Set on the call to MA77_open.
    integer(short) :: npage(2) ! npage for of01
    integer(short) :: nsup ! Number of supervariables, including one for
!           indices that are not used.
    integer(short) :: ntwo ! Number of 2x2 pivots passed by user to analyse
    integer(short) :: null ! Number of null rows and cols in matrix
    logical :: pos_def ! .true. iff the problem is known to positive definite.
!           Specified by user on call of MA77_factor.
    integer(long)  :: posfac ! Holds the position in the main real superfile of
!           the start of the reals in the factor.
    integer(long)  :: posint ! Holds the position in the main int. superfile
!           of the start of the integers in the factor
!           (only needed in the indef. case)
    integer(long)  :: rfree ! First free location in the main real superfile.
    integer(long)  :: rtopmx  ! Max. size of main real stack.
    integer(long)  :: rtopdmx ! Max. size of real stack for delayed pivots.
    integer(short) :: scale = 0 ! Set to 1 if ma77_factor called with
!           scale present.
    integer(short) :: status = 0 ! Monitors the progress:
!      0 initially or after a call to MA77_finalise
!      1 after call to MA77_open
!      2 after a successful call to MA77_analyse
!      3 after a successful call to MA77_factor
!        or successful call to MA77_restart
!     -2 after an error condition early in MA77_open (needed by finalise)
!     -1 after another error condition
!              NB To continue after an error requires a call of MA77_finalise
    integer(long) :: used ! amount of storage used when user wishes to work
!           in-core (ie the sum of the sizes of the arrays
!           that are used in place of files)
    integer(short) :: tnode ! No. of nodes in the tree (including any that are
!           discarded but not reused).

!!! allocatable components:

    real(wp), allocatable :: aelt(:) ! Used in MA77_input_reals in element
!            entry case to squeeze out duplicates/out-of-range entries.
!            Allocated by MA77_input_reals.
!            Deallocated by MA77_factorize.

    real(wp), allocatable :: arow(:) ! Used in MA77_input_reals in row
!            entry case to squeeze out duplicates/out-of-range entries.
!            Allocated by MA77_input_reals.
!            Deallocated by MA77_factorize.

    integer(short), allocatable :: clist(:) ! Used in MA77_input_vars to hold
!            variable list for the incoming
!            element/row with duplicates/out-of-range indices removed.
!            Also used in MA77_input_reals to hold this list
!            Allocated by MA77_input_vars.
!            Deallocated by MA77_factorize.

    integer(long), allocatable :: ifile(:)
!            ifile(ie) is the position in the main integer superfile of the
!            list of names of variables associated with node ie.

    integer(short), allocatable :: iptr(:) ! Used in MA77_input_reals to hold
!            positions of columns in compressed element. 
!            Allocated by MA77_input_vars.
!            Deallocated by MA77_factorize.

    integer(short), allocatable :: map(:) ! Used to hold mapping
!            between user-supplied
!            variable list for incoming element/row  and list
!            with duplicates/out-of-range indices removed.
!            Allocated by MA77_MA77_input_vars.
!            Deallocated by MA77_analyse.

    integer(short), allocatable :: new(:)
!            new(is) is the new supervariable constructed from variables of
!            supervariable is. Allocated by MA77_open to have length n+1.
!            Deallocated by MA77_analyse.

    integer(long), allocatable :: rfile(:)
!            rfile(ie) is the position in the main super-real file after the
!            last entry in the list of reals associated with the input
!            element/row ie.
!            Deallocated by MA77_finalise.

    integer(short), allocatable :: roots(:)
!            Holds the roots of the forest.

    integer(short), allocatable :: size(:)
!            size(ie) holds the length of the list of variables associated
!            with node ie. If there are duplicates/out-of-range indices in the
!            lists supplied by the user at leaf nodes, this is indicated
!            by a negative flag and -size(ie) holds no. of variables at leaf ie
!            after duplicates absorbed and out-of-range entries removed.
!            Not altered by factor or solve.

    integer(short), allocatable  :: size_ind(:)
!            size_ind(ie) holds the length of the list that is written to the
!            delayed pivots superfile for node ie (indefinite case only).

    integer(short), allocatable :: splitp(:)
!            splitp(node) best split point for node (number of
!            children of node to be processed before assembly begins)
!            Note : could be held as file data

    integer(short), allocatable :: svar(:)
!            svar(i) holds the supervariable to which variable i belongs.
!            Allocated by MA77_open to have length n.
!            Set during the calls to MA77_input_vars. used by MA77_analyse.
!            note: not altered so that analyse can be repeated without
!            recalling MA77_input_vars we could save storage by writing
!            it to the main integer superfile during factorize and read
!            back at end of factorize so ready to recall analyse?

    integer(short), allocatable :: vars(:)
!            vars(is) holds the number of variables in supervariable is.
!            Note that vars(1) holds number of variables in supervariable 1,
!            that is, the number of indices not used to index a variable.
!            Allocated by MA77_open to have length n+1.
!            Set during the calls to MA77_input_vars. Deallocated by
!            MA77_analyse.

    integer(short), allocatable :: varflag(:)
!            In MA77_input_vars, varflag(is) = ielt if supervariable is has
!            been encountered in element ielt. Allocated by MA77_open
!            to have length n+1. Deallocated by MA77_analyse.

! The following array components are only used if the user wants to work
! in-core; otherwise they are allocated to have length 1.

    integer(long) :: size_imain ! size of keep%imain
    integer(long) :: size_rmain ! size of keep%rmain
    integer(long) :: size_rwork ! size of keep%rwork
    integer(long) :: size_rwdelay ! size of keep%rwdelay

    integer(short),allocatable :: imain(:)
!            main integer array (replaces main integer file)

    real(wp),allocatable :: rmain(:)
!            main real array (replaces main real file)

    real(wp),allocatable :: rwork(:)
!            main real work array (replaces main real work file)

    real(wp),allocatable :: rwdelay(:)
!            real array for delayed pivots (replaces file)

    character,allocatable :: file1(:)
    character,allocatable :: file2(:)
    character,allocatable :: file3(:)
    character,allocatable :: file4(:)
    character,allocatable :: file5(:)

    type (MA77_node),allocatable :: tree(:) ! holds information at
!        non-leaf nodes. Not altered during factor or solve.
! structures for of01
    type (of01_rdata) :: rdata
    type (of01_idata) :: idata

  end type MA77_keep

contains

!****************************************************************************

  subroutine MA77_open_double(n,filename,keep,control,info,nelt,path)
! User must input basic problem data. Data is checked and then
! superfiles are opened using hsl_of01.

    integer(short), intent (in) :: n ! Order of the system. Must be >= 0.
    character (len=*), intent (in) :: filename(4) ! len(filename) < 400
!            filename(1:2) must identify the integer and real superfiles
!            for the factors;
!            filename(3)   must identify the primary real work superfile.
!            filename(4:5) must identify the integer and real superfiles
!            for delayed pivots (indefinite case only).
    type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
    type (MA77_control), intent (in) :: control ! See derived-type declaration
    type (MA77_info), intent (out) :: info      ! See derived-type declaration
    integer(short), optional, intent (in) :: nelt ! If absent entry is by rows;
!          otherwise, must be set to largest integer used to index an element.
    character (len=*), optional, intent (in) :: path(:) ! len(path) < 400.
!          path(:) must be set to hold path names for the direct access
!          files. More than one name may be needed for large problems.
!          The file names are essentially the concatentation of path and
!          filename (full details in of01 spec).
!          if absent, path = '' used.

! Local variables
    integer(short) :: file ! Indicates to error_open which file
!                            (if any) is involved.
    integer(short) :: l    ! Temporary
    integer(long)  :: llong ! Temporary
    integer(long)  :: maxstore ! Copy of control%maxstore
    integer(short) :: nout ! Unit for printing of messages (-1 for no messages)
    integer(short) :: st   ! stat variable for allocate and deallocate
    integer(long)  :: storage(3)

! Possible error returns:
!  -1   Allocation error
!  -3   Call is not first call or does not follow call to finalise.
!  -4   n  <  0
!  -5   Error from of01_open (error in Fortran inquire).
!  -7   Error from of01_open (error in Fortran open).
! -12   Error from of01_open (a file with same name already exists)
! -13   len(filename) too long
! -16   len(path) too long or open unsuccessful for all elements of path
! -18   nelt < 0
! -23   control%buffer_lpage(:) < 1 or > control%file_size
! -26   control%buffer_npage(:) < 1
! -27   control%maxstore out of range

! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      if (present(nelt)) then
        write (control%unit_diagnostics,'(//a,2(/a,i12))') &
     ' Entering MA77_open (element entry) with:', &
     ' n                   Order of matrix                     = ',n, &
     ' nelt                Number of elements                  = ',nelt
      else
        write (control%unit_diagnostics,'(//a,/a,i12)') &
     ' Entering MA77_open (row entry) with:', &
     ' n                   Order of matrix                     = ',n
      end if
      write (control%unit_diagnostics,'(2a)') ' filename(1)     =   ', &
      trim(filename(1))
      write (control%unit_diagnostics,'(2a)') ' filename(2)     =   ', &
      trim(filename(2))
      write (control%unit_diagnostics,'(2a)') ' filename(3)     =   ', &
      trim(filename(3))
      write (control%unit_diagnostics,'(2a)') ' filename(4)     =   ', &
      trim(filename(4))
      if (present(path)) &
      write (control%unit_diagnostics,'(2a)') ' path(1)         =   ', &
      trim(path(1))

      write (control%unit_diagnostics, &
     '(a,4(/a,i12),2(/a,2i8),2(/a,es12.4)/a/3es12.4)') &
     ' control parameters (control%) :', &
     ' print_level         Level of diagnostic printing           = ', &
       control%print_level, &
     ' unit_diagnostics    Unit for diagnostics                   = ', &
       control%unit_diagnostics, &
     ' unit_error          Unit for errors                        = ', &
       control%unit_error, &
     ' unit_warning        Unit for warnings                      = ', &
       control%unit_warning, &
     ' buffer_lpage        Length of pages in integer/real buffer = ', &
       control%buffer_lpage(1:2), &
     ' buffer_npage        Number of pages in integer/real buffer = ', &
       control%buffer_npage(1:2), &
     ' file_size           Target size for each file              = ', &
       real(control%file_size), &
     ' maxstore            Storage for in-core arrays             = ', &
       real(control%maxstore), &
     ' storage(1:3)        Initial sizes for in-core arrays       = ', &
       real(control%storage(1:3))

    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    file = 0

!    write (6,*) ' on open, status = ',keep%status

! Check status parameter
    if (keep%status /= 0) then
      info%flag = -3
      call error_open; return
    end if
    keep%status = -2

    keep%lup = huge(0_long)/8
    if (control%bits == 32) keep%lup = huge(0_short)/8

! Perform tests on control parameters/ input data
    if (control%buffer_lpage(1) < 1 .or. &
        control%buffer_lpage(1) > control%file_size) then
      info%flag = -23
    else if (control%buffer_lpage(2) < 1 .or. &
        control%buffer_lpage(2) > control%file_size) then
      info%flag = -23
    else if (len(filename) > 400) then
      info%flag = -13
    else if (control%buffer_npage(1) < 1 .or. control%buffer_npage(2) < 1) then
      info%flag = -26
    else if (control%maxstore < 0 .or. control%maxstore > keep%lup) then
      info%flag = -27
    else if (n < 0) then
      info%flag = -4
    end if

    if (info%flag < 0) then
      call error_open; return
    end if

    if (present(nelt)) then
      if (nelt < 0) info%flag = -18
    end if
    if (present(path)) then
! The user must supply file names of length <400 and
! path names of length <400
      if (len(path) > 400) info%flag = -16
    end if

    if (info%flag < 0) then
      call error_open; return
    end if

    keep%n = n
    keep%nelt = n
    l = 1 + n*2
    keep%element_input = .false.
    if (present(nelt)) then
      keep%nelt = nelt
      l = 1 + nelt*3
      keep%element_input = .true.
    end if

! Separate deallocate statements are used so that all are deallocated even
! if some are already deallocated.
    deallocate (keep%imain,stat=st)
    deallocate (keep%rmain,stat=st)
    deallocate (keep%rwork,stat=st)
    deallocate (keep%rwdelay,stat=st)
    deallocate (keep%size,stat=st)

! Allocate the arrays needed by MA77_input_vars
    deallocate (keep%varflag,stat=st)
    deallocate (keep%vars,stat=st)
    deallocate (keep%svar,stat=st)
    deallocate (keep%new,stat=st)
    deallocate (keep%rfile,stat=st)
    deallocate (keep%ifile,stat=st)

    allocate (keep%varflag(n+1),keep%vars(n+1), &
      keep%svar(n),keep%new(n+1),keep%ifile(l), &
      keep%rfile(l),keep%size(l),stat=st)
    if (st /= 0) then
      info%flag = -1; info%stat = st
      call error_open; return
    end if

! Initialise the leaf nodes of the tree
    l = keep%nelt+1
    keep%size(1:l)  = 0
    keep%ifile(1:l) = 0

    keep%rfile = 0

! Set keep%svar and keep%vars to represent all variables
! belonging to supervariable 1. Initialize keep%nsup and keep%varflag.
    keep%svar(:) = 1
    keep%vars(1) = n
    keep%varflag(:) = 0
    keep%nsup = 1

! Initialise of01 (we do this even if user hopes to work in-core)
    keep%status    = -1
    keep%file_size = control%file_size
    keep%npage     = control%buffer_npage
    keep%lpage     = control%buffer_lpage

    if (present(path)) then
      call of01_initialize(info%flag,keep%idata,path=path,npage=keep%npage(1),&
           lpage=keep%lpage(1),file_size=keep%file_size,lp=-1)
    else
      call of01_initialize(info%flag,keep%idata,npage=keep%npage(1), &
           lpage=keep%lpage(1),file_size=keep%file_size,lp=-1)
    end if
    if (info%flag == -1) then
      info%stat = keep%idata%stat
      call error_open;  return
    end if

    if (present(path)) then
      call of01_initialize(info%flag,keep%rdata,path=path,npage=keep%npage(2),&
           lpage=keep%lpage(2),file_size=keep%file_size,lp=-1)
    else
      call of01_initialize(info%flag,keep%rdata,npage=keep%npage(2), &
           lpage=keep%lpage(2),file_size=keep%file_size,lp=-1)
    end if
    if (info%flag == -1) then
      info%stat = keep%rdata%stat
      call error_open;  return
    end if

! integer superfile
    file = 1
    call of01_open(filename(1),keep%index(1),info%flag,keep%idata,lp=-1)
    if (info%flag < 0) then
      info%iostat = keep%idata%iostat
      call error_open;  return
    end if

    do file = 2,4
! real superfiles
      call of01_open(filename(file),keep%index(file),info%flag,&
           keep%rdata,lp=-1)
      if (info%flag < 0) then
        info%iostat = keep%rdata%iostat
        call error_open; return
      end if
    end do

! Check control%maxstore and control%storage(1:3) to see if the user
! wants to try in-core working. Note: at this stage, not known if the
! problem is pos. def. or indef.
    keep%used = 0_long
    maxstore = control%maxstore
    if (maxstore > 0_long) then
      storage(1:3) = control%storage(1:3)
! user wants to try in-core working
      if (any(storage(1:3) == 0_long)) then
        llong = max(1_long,maxstore/200_long)
        storage(1) = 2*llong; storage(2) = llong*8; storage(3) = llong
      end if
      if (storage(2) > 0 .and. keep%used + 2*storage(2) <= maxstore) then
        deallocate (keep%rmain,stat=st)
        allocate (keep%rmain(storage(2)),stat=st)
        if (st == 0) then
! array allocated for real data so we will start by working in-core.
! set unit number for real superfile to be negative to indicate this
          keep%index(2) = -keep%index(2)
          keep%used = keep%used + 2*storage(2)
          keep%size_rmain = storage(2)
        end if
      end if
      if (storage(1) > 0_long .and. storage(1) <= maxstore) then
        deallocate (keep%imain,stat=st)
        allocate (keep%imain(storage(1)),stat=st)
        if (st == 0) then
! array allocated for integer data so we will start by working in-core.
! set unit number for integer superfile to be negative to indicate this
          keep%index(1) = -keep%index(1)
          keep%used = keep%used + storage(1)
          keep%size_imain = storage(1)
        end if
      end if

      if (storage(3) > 0_long .and. keep%used + 2*storage(3) <= maxstore) then
        deallocate (keep%rwork,stat=st)
        allocate (keep%rwork(storage(3)),stat=st)
        if (st == 0) then
          keep%index(3) = -keep%index(3)
          keep%used = keep%used + 2*storage(3)
          keep%size_rwork = storage(3)
        end if
      end if
    end if

    keep%maxstore = maxstore
    keep%flag = 0
    keep%status = 1
    keep%inelrs = 0
    keep%inelrn = 0
    keep%maxfront = 0
    keep%maxlen = 0
    keep%mvar = 0
    keep%mmvar = 0
    keep%matrix_dup = 0
    keep%matrix_outrange = 0

    keep%ifree = 1
    keep%rfree = 1
    keep%posfac = 1

    info%flag = 0
    info%index(1:4) = abs(keep%index(1:4))
    info%nsup = 1
    info%storage(1:4) = 0_long

    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
      write (control%unit_diagnostics,'(/a,i4/a,4i4)') &
     ' Leaving MA77_open with error flag        info%flag  = ', info%flag, &
     ' File indices                        info%index(1:4) = ',info%index(1:4)

  contains

    subroutine error_open
      integer(long)  :: lenw
      integer(short) :: flag,num_file

      keep%name = 'MA77_open'

      if (info%flag == -5 .or. &
          info%flag == -7 .or. info%flag == -12 .or. info%flag == -17) then
! close any files that had already been opened
        keep%status = -2
        if (file > 1) call of01_close(abs(keep%index(1)),lenw,num_file,flag, &
           keep%idata,lp=-1,lkeep=.false.)
        if (file > 2) call of01_close(abs(keep%index(2)),lenw,num_file,flag, &
           keep%rdata,lp=-1,lkeep=.false.)
        if (file > 3) call of01_close(abs(keep%index(3)),lenw,num_file,flag, &
           keep%rdata,lp=-1,lkeep=.false.)
        call of01_end(flag,keep%idata,-1)
        call of01_end(flag,keep%rdata,-1)
      end if

! If necessary, reset hsl_of01 flag value to correct hsl_ma77 value.
! Then print error message.
      if (info%flag == -17) info%flag = -16
      call MA77_print_iflag(keep,nout,info%flag)
      if (info%flag == -12 .and. nout >= 0) &
        write (nout,*) 'filename = ', trim(filename(file))

    end subroutine error_open

  end subroutine MA77_open_double

!****************************************************************************

  subroutine MA77_input_vars_double(index,nvar,list,keep,control,info)

! This subroutine is called by the user to specify the variables
! belonging to an element/row (reverse communication interface).
! Lists of supervariables are generated.
! Note: supervariable 1 holds the variables that are not
! used (it may have no variables belonging to it).
! It is possible to have nsup = n+1.
! Duplicates and/or out of range entries are allowed.

    integer(short), intent (in) :: index ! Index of incoming element/row.
    integer(short), intent (in) :: nvar ! Number of variables
!            in the incoming element/row. Must be >= 0.
    integer(short), intent (in) :: list(nvar) ! Names of the variables in the
!            incoming element/row.
    type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
    type (MA77_control), intent (in) :: control ! See derived-type declaration
    type (MA77_info), intent (inout) :: info    ! See derived-type declaration

! Local scalars
    integer(short) :: flag     ! local error flag
    integer(short) :: i        ! Temporary
    integer(short) :: is       ! Supervariable
    integer(short) :: js       ! New supervariable constructed from
!           variables of is
    integer(short) :: jvar     ! Index for a variable
    integer(short) :: k        ! Temporary
    integer(short) :: l        ! Temporary
    integer(short) :: n1       ! Holds no. of out-of-range entries
    integer(short) :: n2       ! Holds no. of duplicated entries
    integer(short) :: nnvar    ! Set to nvar - n1 - n2 (no. of entries
!           in variable list after duplicates/out-of-range entries removed).
    integer(short) :: nout     ! Unit for printing of error messages
    integer(short) :: nout1    ! Unit for printing of warning messages
    integer(short) :: nsup     ! No. of supervariables
    integer(short) :: st       ! Allocation error
    integer(long) :: loc       ! Location in main integer superfile
    integer(short) :: ldiag

! Possible error returns (info%flag):
! -1   Allocation error
! -3   Error in sequence of calls to routines in package. Immediate return.
! -5   Error from of01_write (error in Fortran inquire)
! -7   Error from of01_write (error in Fortran open)
! -9   MA77_input_vars has already been called for this element/row.
! -15  Error from of01_write (error in Fortran write)
! -31  All entries in element/row are out-of-range
! -33  nvar < 0. No action taken.

    if (keep%n == 0) return
! Perform appropriate printing
    if (control%print_level >= 2 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a,2(/a,i12))') &
     ' Entering MA77_input_vars with:', &
     ' index               Index of element/row                = ',index, &
     ' nvar                Number of variables                 = ', nvar
      i = min(nvar,10)
      write (control%unit_diagnostics,'(a,2(/5i12))') &
     ' Input index list: ',list(1:i)
      if (i < nvar) write (control%unit_diagnostics,'(a)') '  . . . . . .'
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    nout1 = control%unit_warning
    if (control%print_level < 0) nout1 = -1

    info%flag = keep%flag
    keep%name = 'MA77_input_vars'

    if (keep%status /= 1) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      if (keep%status == 0) then
! MA77_open was not called
         keep%status = -2
      else
         keep%status = -1
      end if
      return
! Check input data
    else if (index < 1 .or. index > keep%nelt) then
      return
    else if (nvar < 0) then
      info%flag = -33;  keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      return
    else if (keep%size(index) /= 0) then
      info%flag = -9;  keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      return
    end if

! Treat the special case nvar = 0 (nothing to do)
    if (nvar == 0) then
      keep%size(index) = 0
      keep%ifile(index) = keep%ifree
      keep%rfile(index) = -1
      return
    end if

    keep%mvar = max(keep%mvar,nvar)

    st = 0
    if (.not. allocated(keep%clist)) then
       allocate (keep%clist(keep%mvar+1),keep%map(0:keep%mvar),stat=st)
    else if (keep%mvar+1 > size(keep%clist)) then
       deallocate (keep%clist,stat=st)
       deallocate (keep%map,stat=st)
       allocate (keep%clist(keep%mvar+1),keep%map(0:keep%mvar),stat=st)
    end if
    if (st /= 0) then
      info%stat = st
      info%flag = -1; keep%status = -1
      if (nout >= 0) &
        call MA77_print_iflag(keep,nout,info%flag,ie=index,st=info%stat)
      deallocate (keep%clist,keep%map,stat=st)
      return
    end if
    keep%clist(1:nvar) = list(1:nvar)

    n1 = 0;   n2 = 0
! Loop over variables in current element/row, checking for out of
! range entries and/or duplicates (both of which result in a warning).
! Use the sign of keep%svar to flag whether a variable has appeared
! in the current element/row.
    l = 0
    ldiag = 1
    do i = 1, nvar
      jvar = keep%clist(i)
      if (jvar < 1 .or. jvar > keep%n) then
! Out of range
        keep%map(i) = 0
        n1 = n1 + 1
        if (control%print_level > 1 .and. nout1 >= 0 &
            .and. keep%matrix_outrange+n1 <= 10) &
          write (nout1,'(a,i8,a,i5,a,i8)') &
        ' Warning: In element/row ', index, '   Entry ', i,&
        '   in list has value ',jvar
        if (info%flag == 0 .or. info%flag == 2) then
          info%flag = info%flag + 1
          call MA77_print_iflag(keep,nout1,info%flag,ie=index)
        end if
      else if (keep%svar(jvar) > 0) then
        keep%svar(jvar) = -keep%svar(jvar)
        l = l + 1
        keep%clist(l) = jvar
        keep%map(i) = l
      else
! We have a duplicate
        n2 = n2 + 1
        if (info%flag == 0 .or. info%flag == 1) then
          info%flag = info%flag + 2
          call MA77_print_iflag(keep,nout1,info%flag,ie=index)
        end if
        do k = 1,l
          if (jvar == keep%clist(k)) exit
        end do
        keep%map(i) = k
        if (control%print_level > 1 .and. nout1 >= 0 &
            .and. keep%matrix_dup+n2 <= 10) then
          write (nout1,'(a,i8/a,i8,a,i8,a,i8)') &
        ' Warning: In element/row ', index, ' Variable ', jvar,&
        '   is repeated in positions ',i,' and ',k
        end if
      end if
      if (index == jvar) ldiag = 0
    end do
    if (keep%element_input) then
      ldiag = 0
    else if (ldiag == 1) then
! have to add diagonal to list (row entry only)
      keep%clist(l+1) = index
      keep%svar(index) = -keep%svar(index)
    end if
!    write (6,*) 'index,ldiag,clist',index,ldiag,keep%clist(1:l+ldiag)

! The number of variables in compressed list is now l
    nnvar = l
    keep%size(index) = nnvar
! use a negative sign to indicate out-of-range/duplicates have been found
    if (n1 > 0 .or. n2 > 0) then
      keep%size(index) = -nnvar
! Update the count of duplicates and out-of-range entries and issue warning
      keep%matrix_outrange = keep%matrix_outrange + n1
      keep%matrix_dup = keep%matrix_dup + n2
! Error if ALL entries out-of-range
      if (n1 == nvar) then
        info%flag = -31;  keep%status = -1
        if (nout >= 0) &
          call MA77_print_iflag(keep,nout,info%flag,ie=index)
        deallocate (keep%clist,keep%map,stat=st)
        return
      end if
    end if

! Loop over variables in compressed list, decrementing the counts of variables
! in supervariables and restoring signs of keep%svar
      do i = 1, nnvar+ldiag
        jvar = keep%clist(i)
        is = -keep%svar(jvar)
        keep%svar(jvar) = is
!       write (6,*) 'i,jvar,is,vars',i,jvar,is,keep%vars(is)
        keep%vars(is) = keep%vars(is) - 1
      end do

! Loop over variables, incrementing the count and resetting keep%vars.
      nsup = keep%nsup
      do i = 1, nnvar+ldiag
        jvar = keep%clist(i)
        is = keep%svar(jvar)
!         write (6,*) 'i,jvar,is,varflag,vars',&
!          i,jvar,is,keep%varflag(is),keep%vars(is)
        if (keep%varflag(is) /= index) then
! First occurrence of supervariable is for current element/row
          keep%varflag(is) = index
          if (keep%vars(is) > 0 .or. is == 1) then
! Establish new supervariable
            nsup = nsup + 1
            keep%vars(nsup) = 1
            keep%varflag(nsup) = index
! new(is) is the new supervariable constructed from
! variables of supervariable is.
            keep%new(is) = nsup
            keep%svar(jvar) = nsup
          else
! No new supervariable needed
            keep%vars(is) = 1
            keep%new(is) = is
          end if
        else
! Subsequent occurrence of supervariable is for current element/row
          js = keep%new(is)
          keep%vars(js) = keep%vars(js) + 1
          keep%svar(jvar) = js
        end if
      end do
      keep%nsup = nsup
!      write (6,*) 'svar',keep%svar(1:keep%n)

! Store compressed variable list in main integer superfile
    loc = keep%ifree
!   write (6,*) ' Stored list: ',keep%clist(1:nnvar)
    call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
         loc,nnvar,keep%clist,flag,keep%idata,-1,keep%maxstore,keep%used)

    if (flag < 0) then
      info%iostat = keep%idata%iostat
      info%flag = flag;  keep%status = -1
      if (nout >= 0) &
        call MA77_print_iflag(keep,nout,info%flag,ie=index,ios=info%iostat)
      deallocate (keep%clist,keep%map,stat=st)
      return
    end if

    if (n1 > 0 .or. n2 > 0) then
! If duplicates/out-of-range, also store nvar plus the mapping from the
! user list into the compressed list
      keep%map(0) = nvar
      loc = keep%ifree + nnvar
!     write (6,*) 'nvar,map',keep%map(0:nvar)
      call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
           loc,nvar+1,keep%map(0:nvar),flag,keep%idata,-1, &
           keep%maxstore,keep%used)

      if (flag < 0) then
        info%iostat = keep%idata%iostat
        info%flag = flag;  keep%status = -1
        if (nout >= 0) &
          call MA77_print_iflag(keep,nout,info%flag,ie=index,ios=info%iostat)
        deallocate (keep%clist,keep%map,stat=st)
        return
      end if
      if (keep%element_input) then
        st = 0
        if (.not. allocated(keep%iptr)) then
          allocate (keep%iptr(keep%mvar),stat=st)
        else if (keep%mvar > size(keep%iptr)) then
          deallocate (keep%iptr,stat=st)
          allocate (keep%iptr(keep%mvar),stat=st)
        end if
        if (st /= 0) then
          info%stat = st
          info%flag = -1; keep%status = -1
          if (nout >= 0) &
            call MA77_print_iflag(keep,nout,info%flag,ie=index,st=info%stat)
          deallocate (keep%clist,keep%map,keep%iptr,stat=st)
          return
        end if
      end if
    end if

! Store position of element/row variable list in the integer superfile
! and set flag to indicate reals not yet stored.
    keep%ifile(index) = keep%ifree
    keep%rfile(index) = -1

! Move pointer for next free location in the integer superfile
    keep%ifree = keep%ifree + nnvar
    if (n1 > 0 .or. n2 > 0) keep%ifree = keep%ifree + nvar + 1
    keep%mmvar = max(keep%mmvar,nnvar)
! In row case, we will allow extra space (in case the diagonal
! has not been included by user as will have to be added in
! when constructing the tree).
    if (.not. keep%element_input) keep%mmvar = max(keep%mmvar,nnvar+1)

    info%unused = keep%vars(1)
    info%nsup   = keep%nsup
    keep%ltree  = min(keep%nsup,2*keep%nelt)

    info%matrix_dup = keep%matrix_dup
    info%matrix_outrange = keep%matrix_outrange

    keep%flag = info%flag

    if (control%print_level >= 2 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a,3(/a,i12))') &
     ' Leaving MA77_input_vars with:', &
     ' Error flag                      info%flag             = ', info%flag, &
     ' Number of duplicated indices    info%matrix_dup       = ',  &
       info%matrix_dup, &
     ' Number of out-of-range indices  info%matrix_outrange  = ', &
       info%matrix_outrange
      i = min(nnvar,10)
      write (control%unit_diagnostics,'(a,2(/5i12))') &
     ' Stored index list: ',keep%clist(1:i)
      if (i < nnvar) write (control%unit_diagnostics,'(a)') '  . . . . . .'
    end if

  end subroutine MA77_input_vars_double

!****************************************************************************

  subroutine MA77_analyse_double(order,keep,control,info)

! Analyse phase.
! The user inputs the pivot order and, using the information
! collected on the calls to MA77_input_vars, the assembly tree is constructed.

! For details of keep, control, info : see derived type descriptions
    type (MA77_keep), intent (inout) :: keep
    integer(short), intent (inout), dimension(keep%n) :: order
!        If i is used to index a variable, |order(i)| must
!        hold its position in the pivot sequence. If a 1x1 pivot i is required,
!        the user must set order(i)>0. If a 2x2 pivot involving variables
!        i and j is required, the user must set
!        order(i)<0, order(j)<0 and |order(j)| = |order(i)|+1.
!        If i is not used to index a variable,
!        order(i) may have any value and this is replaced by zero.
!        On exit, holds the pivot order to be used by MA77_factor.
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info

! Local arrays
    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree
    integer(short), allocatable :: first(:)
!            first(is) is the first element associated with supervariable is.
    integer(short), allocatable :: child(:)
!            List of children of a node in the tree. Copied to tree(ie)%child
    integer(short), allocatable :: iwork(:)
!            Work array used when sorting.
    integer(short), allocatable :: next(:)
!            If ie in use, next(ie) is the next element associated with the
!            same supervariable as element ie or next root if a root.
!            If ie not in use, next(ie) is next node in a list of tree nodes
!            starting from the node with name free or -1 for the final node.
    integer(short), allocatable :: perm(:)
!            perm(i) holds the variable that is i-th in the pivot sequence.
!            Also used for checking user-supplied permutation.
    integer(short), allocatable :: varlist(:)
!            List of variables, as read using of01_read.
    integer(short), allocatable :: wlist(:)
!            Work array used when sorting.
    integer(short), allocatable :: index(:) ! used in ordering nodes
    integer(short), allocatable :: list(:)
!            Array for holding lists when constructing nodes in the tree.
!            This list is written to main integer superfile. Length n+1.
    integer(short), allocatable :: nels(:)
!            nels(i) is the number of elements/rows involving variable i.
    integer(long), allocatable :: child_size(:)
!            child_size(ic) is set to the size of the contribution
!            block for child ic.
    integer(long), allocatable :: m0(:), m1(:)
!            allocated to have length (mx_nchild+1) and used in computing
!            memory needed if the assembly is done after the jth child.
    integer(long), allocatable :: mem(:),child_mem(:)

! Local scalars.
    integer(short) :: cblock ! no. of variables in contribution
!          block for cnode
    integer(short) :: ccnode  ! child node
    integer(short) :: depth ! used to compute depth of tree
    integer(short) :: first_free ! first free node
    integer(short) :: first_root ! First root node
    integer(short) :: flag ! local error flag
    integer(short) :: i ! temporary variable
    integer(short) :: ie ! new generated element
    integer(short) :: ielt ! an element
    integer(long)  :: ifree ! first free location in integer superfile
    integer(short) :: ifs ! position in list of first non fully
                      ! summed variable
    integer(short) :: ii ! temporary variable
    integer(short) :: iswap
    integer(short) :: ir ! do loop variable
    integer(short) :: is ! supervariable that ivar belongs to
    integer(short) :: ivar ! variable to be eliminated
    integer(short) :: j ! temporary variable
    integer(short) :: jcnode  ! child node
    integer(short) :: jj ! temporary variable
    integer(short) :: jdum ! do loop variable
    integer(short) :: jelt ! an element
    integer(short) :: jvar ! a variable in an element variable list
    integer(short) :: jjvar ! a variable
    integer(short) :: js ! supervariable
    integer(short) :: k ! temporary variable
    integer(short) :: ks ! supervariable
    integer(short) :: k1 ! temporary variable
    integer(short) :: kk ! temporary variable
    integer(short) :: kpiv ! number of variables that have been ordered so far
                    ! (incremented within elim_order)
    integer(short) :: kvar ! a variable that is the pair of ivar in 2x2 pivot
    integer(short) :: l ! temporary variable
    integer(short) :: llist ! length of array list (=n+1)
    integer(long)  :: lnj ! used to compute flop count
    integer(long)  :: loc ! position in integer superfile
!    integer(long)  :: mem ! memory needed for the subtree rooted at root
    integer(short) :: mx_nchild ! largest number of children a node has
    integer(short) :: n ! order of linear system
    integer(short) :: nc !  number of children of a node
    integer(short) :: nchild ! number of children of a node ie
    integer(short) :: nelim ! number of variables eliminated at a node
    integer(short) :: nelim_i ! no. of variables eliminated at node ielt
    integer(short) :: nelt ! copy of keep%nelt
                       ! (tree(ielt)%nelim)
    integer(short) :: nemin ! min. number of eliminations (control%nemin)
    integer(short) :: nfs ! the number of non fully summed variables at a node
    integer(short) :: node ! node in tree
    integer(short) :: nout ! unit for error messages
    integer(short) :: nroot ! no. of roots (ie no. of components)
    integer(long)  :: nschur ! used to compute entries in factor
    integer(short) :: pnode ! parent node in tree
    integer(short) :: nvar ! number of variables in element ielt
    integer(short) :: root ! a root of tree
    integer(short) :: st ! stat parameter
    logical :: belong !

! Possible error returns:
! -1   Allocation error
! -3   Error in sequence of calls to routines in package
! -5   Error from of01_read (error in Fortran inquire)
! -6   Error from of01_read (error in Fortran read)
! -7   Error from of01_read (error in Fortran open)
! -16  Error -17 returned from of01_read
! -21  Error in order. Immediate return.

    if (keep%n == 0) then
      keep%status = 2; return
    end if

! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
     ' Entering MA77_analyse with:'
      write (control%unit_diagnostics,'(a,5(/a,i12))') &
     ' control parameters (control%) :', &
     ' print_level         Level of diagnostic printing        = ', &
       control%print_level, &
     ' unit_diagnostics    Unit for diagnostics                = ', &
       control%unit_diagnostics, &
     ' unit_error          Unit for errors                     = ', &
       control%unit_error, &
     ' unit_warning        Unit for warnings                   = ', &
       control%unit_warning, &
     ' nemin               Node amalgamation parameter         = ', &
       control%nemin
     if (control%print_level >= 2) then
! Print out pivot order.
       i = min(10,keep%n)
       write (control%unit_diagnostics,'(a,2(/5i12))')  &
     ' User-supplied elimination order :', order(1:i)
       if (i < keep%n) write (control%unit_diagnostics,'(a)') '  . . . . . .'
      end if
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1

    info%flag = keep%flag
    keep%name = 'MA77_analyse'
    ifree = keep%ifree
    info%nfactor = 0_long
    info%nflops =  0_long
    info%ntwo = 0
    keep%ntwo = 0
    keep%maxelim = 0

    info%nio_read(1) = keep%idata%nio_read
    info%nio_read(2) = keep%rdata%nio_read
    info%nio_write(1) = keep%idata%nio_write
    info%nio_write(2) = keep%rdata%nio_write
    info%nwd_read(1) = keep%idata%nwd_read
    info%nwd_read(2) = keep%rdata%nwd_read
    info%nwd_write(1) = keep%idata%nwd_write
    info%nwd_write(2) = keep%rdata%nwd_write

! Deallocate arrays used by MA77_input_vars that are no longer needed
    deallocate (keep%new,stat=st)
    deallocate (keep%vars,stat=st)
    deallocate (keep%map,stat=st)

! Check status parameter (remember that input of integer and real
! data can be interleaved)
    if (keep%status /= 1) then
      info%flag = -3
      call error_anal;  go to 500
    end if

    n = keep%n
    nemin = control%nemin
! Check nemin (a node is merged with its parent if both involve
! fewer than nemin eliminations). If out of range, use the default
    if (nemin < 1) nemin = nemin_default

! Check the user-supplied array order and set the inverse in perm.
! Also add up number of variables that are not used (null rows)
    deallocate (perm,stat=st)
    allocate (perm(n),stat=st)
    if (st /= 0) go to 490
    perm(:) = 0
    k1 = 0
    keep%null = 0
    do i = 1, n
       if (keep%svar(i) == 1) then
! Variable i is not used
          order(i) = 0
          keep%null = keep%null + 1
       else
          j = order(i)
          jj = abs(j)
          if (jj < 1 .or. jj > n) exit
          if (perm(jj) /= 0) exit ! Duplicate found
          if (j > 0) then
            perm(jj) = i
          else
            k1 = k1 + 1
            perm(jj) = -i
          end if
       end if
    end do
    if (i-1 /= n) then
      info%flag = -21
      call error_anal;  go to 500
    end if

! If 2x2 pivots entered, check that pairs are adjacent in perm (error if not)
! and set the second entry in the pair to have positive flag
    if (k1 > 0) then
      l = 1
      do
        if (l > n) exit
        k = perm(l)
        if (k < 0) then
          kk = perm(l+1)
          if (kk >= 0) then
! 2x2 pivot not adjacent
            info%flag = -21
            call error_anal;  go to 500
           end if
           keep%ntwo = keep%ntwo + 1
           perm(l+1) = -kk ! (so the second in the pair is positive)
           l = l + 2
         else
           l = l + 1
         end if
       end do
     end if
!     write (6,*) 'order',order(1:n)
!     write (6,*) 'perm ',perm(1:n)

    keep%l1 = keep%nelt + 1
    keep%l2 = keep%nelt + keep%ltree

! Set up space for the tree (no space needed for the leaf nodes)
    deallocate (keep%tree,stat=st)
    allocate (keep%tree(keep%l1:keep%l2),stat=st)
    if (st /= 0) go to 490

! Allocate varlist (for holding lists of variables)
    allocate (varlist(keep%mmvar),stat=st)
    if (st /= 0) go to 490
! Allocate nels to hold number of elements containing each variable. Initialise
    allocate (nels(n),stat=st)
    if (st /= 0) go to 490
    nels(:) = 0

! Set up an empty list of unused nodes and root nodes (links are in next(:))
    first_free = -1
    first_root = -1
    nroot = 0
    keep%tnode = keep%nelt
    keep%maxlen = keep%mmvar
    keep%mx_ifree = 0_long

! Allocate arrays first and next for holding linked lists of elements
! (link according to the supervariable
! in its list that appears first in the pivot sequence
    deallocate (first,stat=st)
    deallocate (next,stat=st)
    allocate (first(keep%nsup),next(keep%l2),stat=st)
    if (st /= 0) go to 490
    first(:) = 0

! keep%varflag has already been allocated by ma77_open. initialise.
    keep%varflag(1:n) = 0

! The array child will hold a temporary list of the children
! of the node under construction.
    llist = n + 1 ! extra space means always one space between variables
                  ! that are fully summed and those that are not
    allocate (child(1),list(llist),stat=st)
    if (st /= 0) go to 490

! Call internal subroutine to construct the tree
    if (keep%ntwo == 0) then
! no 2x2 pivots
      call construct_tree1
    else
      call construct_tree2
    end if
    if (st /= 0) go to 490
    if (info%flag < 0) then
      call error_anal; go to 500
    end if

! Store first free location in main integer superfile where we can write
! lists of integers (will only need to write integer lists in indef. case)
    keep%posint = keep%ifree

! Deallocate arrays we are done with.
    deallocate (child,stat=st)
    deallocate (list,stat=st)
    deallocate (first,stat=st)
    deallocate (keep%varflag,stat=st)

! Store the roots that were found during the tree construction
    deallocate (keep%roots,stat=st)
    allocate (keep%roots(nroot),stat=st)
    if (st /= 0) go to 490
    i = first_root
    j = 0
    do
      if (i <= 0) exit
      j = j + 1
      keep%roots(j) = i
      i = next(i)
    end do
    deallocate (next,stat=st)

!**************************************
! Find the split point for each non-leaf node in the tree
    deallocate (keep%splitp,stat=st)
    deallocate (iwork,stat=st)
    deallocate (child_size,stat=st)
    deallocate (m0,stat=st)
    deallocate (m1,stat=st)
    deallocate (index,stat=st)
    deallocate (mem,stat=st)
    deallocate (child_mem,stat=st)

    nelt = keep%nelt
    allocate (keep%splitp(nelt+1:keep%tnode),child_size(1:mx_nchild),&
       m0(0:mx_nchild),m1(0:mx_nchild),mem(nelt+1:keep%tnode), &
       child_mem(1:mx_nchild),index(1:mx_nchild),iwork(1:mx_nchild),stat=st)
     if (st /= 0) goto 490

     keep%splitp(nelt+1:keep%tnode) = 0
     mem = 0

! First have to ensure, for each node, the list of children is
! organised so that any children that are leaf nodes are at the
! end of the list. Loop over non-leaf nodes.
! In this loop, we also compute info%nfactor and info%nflops
    do node = nelt+1,keep%tnode
      if (.not. allocated(keep%tree(node)%child) ) cycle

      nvar = abs(keep%size(node))
      nelim = keep%tree(node)%nelim

      nschur = nvar - nelim
      nschur = nvar + nschur + 1_long
      lnj = nelim
      info%nfactor = info%nfactor + (lnj*nschur)/2_long
      do j = 0,nelim-1
        lnj = nvar - j
        info%nflops = info%nflops + lnj*lnj
      end do

      nchild = size(keep%tree(node)%child)

! Loop through the children to find the number of non-leaf children
! and swap them to be ahead of leaf children.
      nc = 0
      k1 = nchild
out:  do j = 1, nchild
        jcnode = keep%tree(node)%child(j)
        if (jcnode <= nelt) then
! Leaf node. Look for swap
          do k = k1,j+1,-1
            ccnode = keep%tree(node)%child(k)
            if (ccnode > nelt) then
! non-leaf so swap with jcnode
              nc = nc + 1
              keep%tree(node)%child(j) = ccnode
              keep%tree(node)%child(k) = jcnode
              k1 = k - 1
              cycle out
            end if
          end do
! no swap available so rest of children are leaf nodes
          exit out
        else
          nc = nc + 1
        end if
      end do out

    end do

! Allocate arrays for depth first search of tree
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    allocate (count(0:keep%tnode),cnode(0:keep%tnode),stat=st)
    if (st /= 0) goto 490

    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               info%maxdepth = max(depth,info%maxdepth)
               count(depth) = 0
               cycle
            else
               info%maxdepth = max(depth+1,info%maxdepth)
            end if
         end if

! Perform work at node. 
        call order_child(node)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree

     if (info%flag < 0) then
       call error_anal; go to 500
     end if
   end do
   keep%maxdepth = max(1,info%maxdepth)

!    do ir = 1, nroot
!      root = keep%roots(ir)
!      if (control%print_level == 2) then
!       j = 1; call write_tree(root,0)
!      end if
!    end do

! deallocate arrays we have finished with.
    deallocate (iwork,stat=st)
    deallocate (child_size,stat=st)
    deallocate (m0,stat=st)
    deallocate (m1,stat=st)
    deallocate (index,stat=st)
    deallocate (mem,stat=st)
    deallocate (child_mem,stat=st)

!**************************************

! For each variable i, set order(i) to hold the position
! that i is eliminated at.
    if (size(varlist) < keep%maxlen) then
      deallocate (varlist,stat=st)
      allocate (varlist(1:keep%maxlen),stat=st)
      if (st /= 0) go to 490
    end if

    kpiv = 0
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node.
        call elim_order(order,keep%tree,node,kpiv) 
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) then
        call error_anal; go to 500
      end if
    end do
!   write (6,*) 'after ', order(1:n)


! We now rewrite the index lists for non-leaf nodes in elimination order
    allocate (iwork(keep%maxlen),stat=st)
    if (st /= 0) go to 490
    deallocate (wlist,stat=st)
    allocate (wlist(keep%maxlen),stat=st)
    if (st /= 0) go to 490

    do 270 ielt = nelt+1, keep%tnode
      nvar = keep%size(ielt)
      if (nvar <= 0) cycle
! Cycle if nothing needs sorting
      nelim = keep%tree(ielt)%nelim
      nvar = nvar - nelim
      if (nvar <= 1) cycle
! Read from the main integer superfile the variables that are not eliminated.
      loc = keep%ifile(ielt) + nelim
      call MA77_read_integer(keep%index(1),keep%imain, &
           loc,nvar,varlist,flag,keep%idata,-1)
      if (flag < 0) then
        info%flag = flag
        call error_anal; go to 500
      end if
! Sort variables into ascending order according to when they are eliminated.
      do i = 1, nvar
        jvar = abs(varlist(i))
        k = abs(order(jvar))
        iwork(i) = k
      end do

      call kb07ai(iwork,nvar,wlist)

! Entries of varlist now have to be ordered in the same way.
      do i = 1, nvar
        j = wlist(i)
        iwork(i) = varlist(j)
      end do
! Overwrite the original variable list with the permuted list
      call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
           loc,nvar,iwork,flag,keep%idata,-1,keep%maxstore,keep%used)
      if (flag < 0) then
        info%flag = flag
        call error_anal; go to 500
      end if

270 continue

! Set keep%status to 2 to indicate successful analyse is complete
    keep%status = 2
    info%index(1:4) = keep%index(1:4)
    info%maxfront = keep%maxfront

    info%nio_read(1) = keep%idata%nio_read - info%nio_read(1)
    info%nio_read(2) = keep%rdata%nio_read - info%nio_read(2)
    info%nio_write(1) = keep%idata%nio_write - info%nio_write(1)
    info%nio_write(2) = keep%rdata%nio_write - info%nio_write(2)
    info%nwd_read(1) = keep%idata%nwd_read - info%nwd_read(1)
    info%nwd_read(2) = keep%rdata%nwd_read - info%nwd_read(2)
    info%nwd_write(1) = keep%idata%nwd_write - info%nwd_write(1)
    info%nwd_write(2) = keep%rdata%nwd_write - info%nwd_write(2)

    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a)') &
     ' Leaving MA77_analyse with:'
      write (control%unit_diagnostics, &
     '(a,6(/a,i12),2(/a,es12.4),4(/a,es12.4))') &
     ' information parameters (info%) :', &
     ' flag              Error flag                                   = ', &
       info%flag,     &
     ' maxfront          Forecast maximum frontsize                   = ', &
       info%maxfront, &
     ' ntwo              Number of 2x2 pivots supplied by user        = ', &
       info%ntwo,     &
     ' nsup              Number of supervariables                     = ', &
       info%nsup,     &
     ' tree_nodes        Number of non-leaf nodes in tree             = ', &
       info%tree_nodes,     &
     ' maxdepth          Maximum depth of the tree                    = ', &
       info%maxdepth,     &
     ' nfactor           Forecast number of entries in L              = ', &
       real(info%nfactor), &
     ' nflops            Forecast number of flops                     = ', &
       real(info%nflops), &
 ' nio_read(1)       No. integer records read from disk by OF01_read   = ',&
       real(info%nio_read(1)), &
 ' nio_write(1)      No. integer records written to disk by OF01_write = ',&
       real(info%nio_write(1)), &
     ' nwd_read(1)       No. integer scalars read by OF01_read        = ', &
       real(info%nwd_read(1)), &
     ' nwd_write(1)      No. integer scalars written by OF01_write    = ', &
       real(info%nwd_write(1))

    end if
    go to 500

490 info%flag = -1
    call error_anal

500 continue
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (child,stat=st)
    deallocate (child_size,stat=st)
    deallocate (first,stat=st)
    deallocate (index,stat=st)
    deallocate (iwork,stat=st)
    deallocate (list,stat=st)
    deallocate (nels,stat=st)
    deallocate (m0,stat=st)
    deallocate (m1,stat=st)
    deallocate (mem,stat=st)
    deallocate (child_mem,stat=st)
    deallocate (perm,stat=st)
    deallocate (varlist,stat=st)
    deallocate (wlist,stat=st)
    deallocate (keep%varflag,stat=st)

  contains
!**********************************************
    recursive subroutine write_tree(node,pnode)
      integer(short) :: node,pnode
!!! Write subtree rooted at root whose parent is pnode
!!! The depth is held in j
      integer(short) :: cnode ! child of node
      integer(short) :: i     ! do index
      integer(short) :: nc    ! number of children
!!
       nvar = abs(keep%size(node))
       if (j<=5) &
       write(*,'((a,i6,a,i2,a,i6,a,i4,a,i4))')&
               'node',node,' depth=',j,' parent=',&
                pnode,' nvar=',nvar,' nelim=',keep%tree(node)%nelim
      loc = keep%ifile(node)
      call MA77_read_integer(keep%index(1),keep%imain,&
           loc,nvar,varlist,flag,keep%idata,-1)
      write(*,*)' list:',varlist(1:nvar)
      if (.not.allocated(keep%tree(node)%child))return
!!
!!! Loop over children of node.
      nc = size(keep%tree(node)%child)
      j = j + 1
      do i = 1, nc
        cnode = keep%tree(node)%child(i)
        if (cnode > nelt) then
           call write_tree(cnode,node)
        end if
      end do
      j = j - 1
!!      nvar = abs(keep%size(ielt))
    end subroutine write_tree
!**********************************************

    subroutine error_anal
! Perform actions after an error detected
      info%iostat = keep%idata%iostat
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat,ios=info%iostat)

      if (info%flag == -21) then
! reset keep%ifree (so that it is possible to recall MA77_analyse)
        keep%ifree = ifree
      else
        keep%status = -1
      end if
    end subroutine error_anal
!**********************************************

    subroutine read_varlist(length)
! Read a list of variables from the main integer superfile
    integer(short) :: length ! length of list to be read
      if (size(varlist) < length) then
        deallocate (varlist,stat=st)
        allocate (varlist(1:length),stat=st)
        if (st /= 0) then
          info%flag = -1
          return
        end if
      end if
      call MA77_read_integer(keep%index(1),keep%imain,loc,length,varlist, &
           flag,keep%idata,-1)
      if (flag < 0) info%flag = flag
    end subroutine read_varlist
!**********************************************

    subroutine reallocate_child(k)
      integer(short),intent(in) :: k
! Reallocate the array child from size nchild to size nchild+nc-1
      if (.not. allocated(iwork)) then
        allocate (iwork(1:nchild),stat=st)
      else if (size(iwork) < nchild) then
        deallocate (iwork,stat=st)
        allocate (iwork(1:nchild),stat=st)
      end if
      if (st /= 0) then
        info%flag = -1
        call error_anal; return
      end if
      iwork(1:nchild) = child(1:nchild)

      deallocate (child,stat=st)
      allocate (child(1:k),stat=st)
      if (st /= 0) then
        info%flag = -1
        call error_anal; return
      end if
! Reset children of ie
      child(1:nchild) = iwork(1:nchild)

    end subroutine reallocate_child
!**********************************************

    subroutine elim_order(order,tree,node,kpiv)

! This subroutine sets the elimination order at node
      integer(short), intent (inout), dimension(keep%n) :: order
      type (MA77_node), intent (in) :: tree(keep%l1:keep%l2)
      integer(short), intent (in) :: node
      integer(short), intent (inout) :: kpiv

      integer(short) :: flag  ! error flag
      integer(short) :: i     ! do loop index
      integer(short) :: j     !
      integer(short) :: l     !

! Immediate return for leaf nodes
      if (node <= nelt) return

! Read in list of variables eliminated at node
      loc = keep%ifile(node)
      nelim = tree(node)%nelim
!     write (6,*) 'node,kpiv,nelim',node,kpiv,nelim
      call MA77_read_integer(keep%index(1),keep%imain,loc,nelim,varlist, &
           flag,keep%idata,-1)
      if (flag < 0) then
        info%flag = flag;  return
      end if
! kpiv is the number of variables that have been ordered so far
! (it was initialised to zero before first call to elim_order)
      if (info%ntwo == 0) then
        do i = 1, nelim
          kpiv = kpiv + 1
          j = varlist(i)
          order(j) = kpiv
        end do
      else
        l = 1
        do i = 1, nelim
          if (l > nelim) exit
          kpiv = kpiv + 1
          j = varlist(l)
          if (j > 0) then
            order(j) = kpiv
            l = l + 1
          else
! j is first entry in 2x2 pivot .... flag both entries in order with negative
! sign (this gives the same form as was supplied by the user)
            order(-j) = -kpiv
            kpiv = kpiv + 1
            j = varlist(l+1)
 !!           varlist(l+1) = -j ! set both entries in pivot with negative flag
            order(j) = -kpiv
            l = l + 2
          end if
        end do

        call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain,loc,&
             nelim,varlist,flag,keep%idata,-1,keep%maxstore,keep%used)
        if (flag < 0) then
          info%flag = flag;  return
        end if

      end if

    end subroutine elim_order

!**********************************************

    subroutine order_child(node)
! This subroutine computes the memory needed for the subtree
! rooted at node and finds the split point for processing its
! children (placed in keep%splitp)
     use hsl_kb22_long_integer ! heapsort
     integer(short), intent (in) :: node

     integer(short) :: cnode  ! child node
     integer(short) :: ic
     integer(short) :: nc ! number of non-leaf children of node
     integer(short) :: lc ! last non-leaf candidate
     integer(short) :: inform ! error flag for kb22
     integer(short) :: nchild ! number of children of node
     integer(short) :: nvar ! no. of variables in node
     integer(short) :: p ! best split point
     integer(long)  :: memory ! total memory
     integer(long)  :: minmem ! min. memory
     integer(short) :: store

     if (node <= nelt) return

     nvar = abs(keep%size(node))
     nelim = keep%tree(node)%nelim

     nchild = size(keep%tree(node)%child)
! Deal with the trivial case (only occurs in row entry case) ... must
! ensure splitp is set in case restart is used
! (in which case, keep%splitp(:) is written out so must be defined)
     if (nchild == 0) then
       keep%splitp(node) = 1
       return
     end if

      nc = nchild
! In element case, add up the number of children of node that are leaf nodes
! .... they are the last children (in row case, leaf nodes are original
! rows and they are not part of the tree)
      if (keep%element_input) then
        if (nchild == 0) return
        do i = 1,nchild
          cnode = keep%tree(node)%child(i)
          if (cnode <= keep%nelt) exit
        end do
        nc = i - 1
      end if

! No ordering needed if all children are leaf nodes
     if (nc < 1) then
       keep%splitp(node) = nchild
       return
     end if

     if (nc == 1) then
       keep%splitp(node) = 1
       cnode = keep%tree(node)%child(1)
       mem(node) = mem(cnode)
       return
     end if

! Order the non-leaf children in decreasing order of memory needed.
     do ic = 1, nc
       index(ic) = ic
       cnode = keep%tree(node)%child(ic)
       child_mem(ic) = mem(cnode)
! take copy of child_mem (since we do not want child_mem to be reordered)
       child_size(ic) = child_mem(ic)
     end do
! Form heap (need kb22 because it works with long integers).
     call kb22_build_heap(nc,child_size,inform,index)
     do ic = nc, 1, -1
! Find smallest entry in current heap and put it at end of the list.
       call kb22_get_smallest(ic,child_size,inform,index)
     end do

! Set child_size(ic) to hold storage for the contribution block for child ic.
     do ic = 1, nc
       cnode = keep%tree(node)%child(ic)
! cblock holds no. of variables in the contribution block.
       cblock = abs(keep%size(cnode)) - keep%tree(cnode)%nelim
       child_size(ic) = (cblock*cblock+cblock)/2
     end do

! Find the best split point p (that is, the point at which
! the parent node is to be allocated and assembly done so that min.
! memory is needed). This idea is taken from
! Guermouche and L'Excellent, TOMS 32, 17-32, 2006.
! Loop over the children of node in order of decreasing memory.
! memory is the total memory needed if assembly done after jth child
     store = (nvar*nvar+nvar)/2
     m0(0) = 0
     m1(0) = 0
     minmem = huge(1_long)
     lc = control%p
     if (lc > 2) lc = nc
     do j = 1, lc
! Place child j in right position
       jj = index(j)
       do i = j - 1, 1, -1
         ii = index(i)
!        write (6,*) i,ii,child_mem(ii) - child_size(ii),  &
!                   child_mem(jj) - child_size(jj)
! exit loop if no further swapping needed.
         if (child_mem(ii)-child_size(ii) >= &
             child_mem(jj)-child_size(jj)) exit
         index(i+1) = ii
       end do
       iswap = i + 1
       index(iswap) = jj
       jj = index(j)
       do i = iswap, j
         ii = index(i)
         m1(i) = max(m1(i-1),child_mem(ii)+m0(i-1))
         m0(i) = m0(i-1) + child_size(ii)
       end do
       memory = m1(j)

       if (j < nc) memory = max(memory,store+child_mem(index(j+1)))
! Record the min. total memory ... the best j is the split point p.
! Because of the way we do merges in factorize, we do not set p = 1
! unless nc = 1
       if (memory > minmem .and. control%p==4 .and. j > 2) exit
       minmem = memory
       p = j
     end do

! Store the split point for node and memory needed
     keep%splitp(node) = p
     mem(node) = minmem

! If p = nc, the children are in decreasing order of (child_mem - child_size).
! Otherwise, we must undo the last set of swaps so that the first p
! children only are in decreasing order of (child_mem - child_size).
     if (p /= nc .and. control%p == 4) then
       ii = index(iswap)
       do i = iswap, j - 1
         index(i) = index(i+1)
       end do
       index(j) = ii
     end if

! Do the physical reordering (iwork is used here as a temporary array)
     do j = 1, nc
       jj = index(j)
       iwork(j) = keep%tree(node)%child(jj)
     end do

     keep%tree(node)%child(1:nc) = iwork(1:nc)
!    write(*,*) 'node',node,' has children', iwork(1:nchild)
!     if( nc/=2 .or. p == 1) then
!        write(*,*) 'node',node,' has', nc,  'children', ' p=', p, &
!        ' stack size',minmem,' front size',store
!       do j  = 1, nc
!         cnode = iwork(j)
!         cblock = abs(keep%size(cnode)) - keep%tree(cnode)%nelim
!         cblock = (cblock*cblock+cblock)/2
!         write(*,*)cnode, mem(cnode),cblock, mem(cnode)-cblock
!       end do
!     end if
     

   end subroutine order_child

!**********************************

   subroutine construct_tree1
! Construct tree when all pivots are 1x1

    if (keep%element_input) then
! Link the elements according to the supervariable
! in its list that appears first in the pivot sequence.
! Set nels(ivar) to be number of elements containing variable ivar.
       do ielt = 1, keep%nelt
        nvar = abs(keep%size(ielt))
        if (nvar == 0) cycle
! Read from the main integer superfile
        loc = keep%ifile(ielt)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,varlist, &
             flag,keep%idata,-1)
        if (flag < 0) then
          info%flag = flag; return
        end if
        jj = n + 1
        do i = 1, nvar
          jvar = varlist(i)
          nels(jvar) = nels(jvar) + 1
          if (order(jvar) < jj) then
            k = jvar
            jj = order(k)
          end if
        end do
! ks is supervariable to which k belongs.
        ks = keep%svar(k)
        next(ielt) = first(ks)
        first(ks) = ielt
       end do
    end if
! Note : in row entry case, we count only the generated elements

    mx_nchild = 0 ! will hold largest number of children a node has.
    do 200 jdum = 1, n
      ivar = perm(jdum)
      if (ivar == 0) cycle
! Check whether we have already dealt with this variable.
      if (order(ivar) == 0) cycle
      is = keep%svar(ivar)
!     write (6,'(a,9i4)') 'jdum,ivar,is,keep%nsup', jdum,ivar,is,keep%nsup

      ie = first_free
      if (ie < 0) then
! New node number needed.
        ie = keep%tnode + 1
        keep%tnode = ie
! Initialize node as empty
        keep%size(ie) = 0
        keep%tree(ie)%nelim = 0
        keep%ifile(ie) = 1
        keep%rfile(ie) = 1
      else
! we can reuse a node number that has been freed
        first_free = next(ie)
      end if

! Check array child is large enough to hold children of ie
      nchild = nels(ivar)
      if (size(child) < nchild) then
        deallocate (child,stat=st)
        allocate (child(1:nchild),stat=st)
        if (st /= 0) return
      end if

! Construct the child list for ie.
      ielt = first(is)
!        write (6,*) 'ie,ivar,is,ielt,nels(ivar)',ie,ivar,is,ielt,nels(ivar)
      do i = 1, nels(ivar)
        child(i) = ielt
        ielt = next(ielt)
      end do
!     write (6,*) 'children:',child(1:nchild)

! nelim holds number of variables that can be eliminated at ie.
! The indices of these variables are placed in the array list from the front.
! The indices of the other variables in the front are placed in list from
! the back; ifs holds the position of the first such variable.
      nelim = 0
      ifs = llist + 1

      if (.not.keep%element_input) then
! Row entry.
! Read the variable list for row ivar from the main integer
! superfile into varlist. Check if diagonal present.
        l = 0
        nvar = abs(keep%size(ivar))
        loc = keep%ifile(ivar)
        call read_varlist(nvar)
        if (info%flag < 0) return
        do ii = 1, nvar
          jvar = varlist(ii)
          if (order(jvar) == 0) cycle
! jvar not yet eliminated.
          if (jvar == ivar) l = 1
          js = keep%svar(jvar)
          if (js == is) then
! jvar belongs to the same supervariable as ivar and so
! can be eliminated.
            nelim = nelim + 1
            list(nelim) = jvar
            keep%varflag(jvar) = nelim
          else
! Put jvar at the end of list as not ready for elimination
            ifs = ifs - 1
            list(ifs) = jvar
            keep%varflag(jvar) = ifs
          end if
        end do
        if (l == 0) then
! add in diagonal to list of variables to be eliminated (that is, ivar)
          nelim = nelim + 1
          list(nelim) = ivar
          keep%varflag(ivar) = nelim
        end if
      end if

! Now loop over the children of ie (which are all generated elements
! in the row entry case), reading in the variables
! for each child that are not yet eliminated and adding them into list.
        do i = 1, nchild
          ielt = child(i)
          nelim_i = 0
          if (ielt > keep%nelt) nelim_i = keep%tree(ielt)%nelim
          nvar = abs(keep%size(ielt)) - nelim_i
          loc = keep%ifile(ielt) + nelim_i
          call read_varlist(nvar)
          if (info%flag < 0) return
! varlist holds the variables in element ielt.
            do ii = 1,nvar
              jvar = varlist(ii)
              if (keep%varflag(jvar) == 0) keep%varflag(jvar) = llist+1
            end do
! Merge varlist into list.
            do ii = 1, nvar
              jvar = varlist(ii)
              if (keep%varflag(jvar) == llist+1) then
! First encounter of jvar for element ie.
                js = keep%svar(jvar)
                if (js == is) then
! jvar belongs to the same supervariable as ivar and so can be eliminated
                  nelim = nelim + 1
                  list(nelim) = jvar
                  keep%varflag(jvar) = nelim
                else
! Put jvar at the end of list
                 ifs = ifs - 1
                 list(ifs) = jvar
                 keep%varflag(jvar) = ifs
                end if
              end if

! Decrement nels(jvar)
              nels(jvar) = nels(jvar) - 1
!           write (6,*) 'ie,i,ielt,jvar,nels(jvar)',ie,i,ielt,jvar,nels(jvar)
            end do
        end do

      if (keep%element_input) then
! Element entry.
! Check whether any of the variables in list(ifs:llist) are fully
! summed (that is, nels(jvar) = 0). If they are, move to front of array list.
        k1 = ifs
        do k = k1,llist
          jvar = list(k)
          if (nels(jvar) == 0) then
! jvar can be eliminated so move into front portion of list
            nelim = nelim + 1
            list(nelim) = jvar
            keep%varflag(jvar) = nelim
! make sure we leave no gaps in list
            if (k /= ifs .and. ifs < llist) then
              kvar = list(ifs)
              list(k) = kvar
              keep%varflag(kvar) = k
            end if
            ifs = ifs + 1
          end if
        end do

      end if
!     write (6,*) 'At node= ', ie, 'nelim= ',nelim, &
!   ' list: ',list(1:nelim),list(ifs:llist)
! Set keep%varflag to -ie for variables that are eliminated
! and set > 0 for those that are not eliminated
      do i = 1,nelim
         jvar = list(i)
         order(jvar) = 0
         keep%varflag(jvar) = -ie
      end do
      do i = ifs,llist
         jvar = list(i)
         keep%varflag(jvar) = i
      end do
! nfs is the number of non fully summed variables
      nfs = llist - ifs + 1

      if (nfs > 0) then
! Link new element into list. Loop over the uneliminated variables belonging
! to ie to find the one that occurs earliest in the pivot sequence.
        jj = n + 1
        do i = ifs,llist
          jvar = list(i)
          nels(jvar) = nels(jvar) + 1
!         if (ie >= 40000) write (6,*) 'ie,i,jvar,nels(jvar),jj,order(jvar)',&
!               ie,i,jvar,nels(jvar),jj,order(jvar)
          if (order(jvar) < jj) then
            k = jvar
            jj = order(k)
          end if
        end do
        ks = keep%svar(k)
        next(ie) = first(ks)
        first(ks) = ie
      else
! Node is a root. Link with other roots
        next(ie) = first_root
        first_root = ie
        nroot = nroot + 1
      end if

! reset keep%varflag to zero for variables that are not eliminated
      do i = ifs,llist
         jvar = list(i)
         keep%varflag(jvar) = 0
      end do

! Merge parent ie and child nodes if the list of uneliminated variables
! at the child is the same as the list of variables at the parent
! or if both involve fewer than nemin eliminations.
      k = nchild
      do 80 i = 1, k
        ielt = child(i)
! Cycle if ielt is a leaf node
        if (ielt <= keep%nelt) cycle
! Skip to next child if merging is not needed.
        nelim_i = keep%tree(ielt)%nelim

        if (nelim+nfs > keep%size(ielt)-nelim_i .and. &
           (nelim >= nemin .or. nelim_i >= nemin)) cycle

! Merge. Read in variable list. Only need the first nelim_i entries for ielt.
! Read directly into list
        if (nelim_i /= 0) then
          loc = keep%ifile(ielt)
          call MA77_read_integer(keep%index(1),keep%imain,loc,nelim_i,    &
               list(nelim+1:nelim+nelim_i),flag,keep%idata,-1)
          if (loc+keep%size(ielt) == keep%ifree) &
! store the largest integer storage needed so far
            keep%ifree = loc
          if (flag < 0) then
            info%flag = flag; return
          end if
          nelim = nelim + nelim_i
        end if

! The children of ielt become children of ie.
        nc = size(keep%tree(ielt)%child)
! ie may have more children than before so check size
! of array child. If necessary, allocate larger array child.
        if (nchild+nc-1 > size(child)) then
          call reallocate_child(nchild+nc-1)
          if (info%flag < 0) return
        end if

! Loop over children of ielt, making them children of ie
        if (nc > 0) then
! The first child replaces ielt in the list of children for ie
          l = i
          do j = 1, nc
            jelt = keep%tree(ielt)%child(j)
            child(l) = jelt
            l = nchild + j
          end do
          nchild = nchild + nc - 1
        else
! Flag that the child has gone.
          child(i) = 0
        end if

! Remove node ielt from keep%tree by deallocating its list of children
! and adding it to the list of free nodes.
        deallocate (keep%tree(ielt)%child,stat=st)
        next(ielt) = first_free
        first_free = ielt
        keep%size(ielt) = 0
        keep%tree(ielt)%nelim = 0
        keep%ifile(ielt) = 1
        keep%rfile(ielt) = 1

! Finished with child ielt
80    continue

! Flush out dummy children
      k = nchild
      nchild = 0
      do i = 1, k
        if (child(i) /= 0) then
          nchild = nchild + 1
          child(nchild) = child(i)
        end if
      end do
      mx_nchild = max(mx_nchild,nchild)

! Copy list of children into tree(ie)%child
      allocate (keep%tree(ie)%child(1:nchild),stat=st)
    ! write (6,*) 'node',ie,' has children',child(1:nchild)
      keep%tree(ie)%child(1:nchild) = child(1:nchild)

! Move list of uneliminated variables forward
      k = nelim
      do i = ifs,llist
         k = k + 1
         list(k) = list(i)
      end do

      keep%size(ie) = k
      keep%maxlen = max(keep%maxlen,k)
      keep%maxfront = max(keep%maxfront,k)
!     write (*,*) 'ie,k,maxlen',ie,k,keep%maxlen
      keep%tree(ie)%nelim = nelim
      keep%maxelim = max(nelim,keep%maxelim)

! Store variable list for ie in main integer superfile.
      loc = keep%ifree
!     write(*,*) 'write',k,' integers',' from position',loc
!     write(*,*) 'ie,nelim,list',ie,nelim,list(1:k)
      call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
           loc,k,list,flag,keep%idata,-1,keep%maxstore,keep%used)
      if (flag < 0) then
        info%flag = flag
        return
      end if

! Store position in the integer superfile of this variable list
      keep%ifile(ie) = loc
! Move pointer for next free location in the integer superfile
      keep%ifree = loc + k
      keep%mx_ifree = max(keep%mx_ifree,keep%ifree)

!!! end of main tree loop
200 continue
    info%tree_nodes = keep%tnode - keep%nelt

   end subroutine construct_tree1

!**********************************

   subroutine construct_tree2
! Construct tree when some pivots are 2x2

!      write (6,*) 'order',order(1:n)
!      write (6,*) 'perm ',perm(1:n)
    if (keep%element_input) then
! Link the elements according to the supervariable
! in its list that appears first in the pivot sequence.

! If ivar is 1x1 pivot or is the first entry in a 2x2 pivot,
! set nels(ivar) to be number of elements containing ivar.
! If kvar is the partner of ivar (i.e. kvar is second entry in 2x2 pivot)
! set nels(kvar) to be number of elements containing kvar but NOT ivar

      do ielt = 1, keep%nelt
        nvar = abs(keep%size(ielt))
        if (nvar == 0) cycle
! Read from the main integer superfile
        loc = keep%ifile(ielt)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,varlist, &
             flag,keep%idata,-1)
        if (flag < 0) then
          info%flag = flag
          return
        end if
! Set flag for variables if ielt
        do i = 1, nvar
          jvar = varlist(i)
          keep%varflag(jvar) = ielt
        end do
        jj = n + 1
        do i = 1, nvar
          jvar = varlist(i)
          nels(jvar) = nels(jvar) + 1
          kk = -order(jvar)
          belong = .false.
          if (kk > 0) then
! jvar is part of 2x2 pivot. Check whether it is 1st or 2nd entry
! in the 2x2 pivot.
! If it is the 2nd entry, check to see if its partner kvar is also in ielt
            if (perm(kk) > 0) then
              kvar = -perm(kk-1)
              if (keep%varflag(kvar) == ielt) then
                nels(jvar) = nels(jvar) - 1
                belong = .true.
              end if
            end if
          end if
          if (abs(order(jvar)) < jj) then
            k = jvar
            if (belong) k = kvar
            jj = abs(order(k))
          end if
        end do
! ks is supervariable to which k belongs.
        ks = keep%svar(k)
        next(ielt) = first(ks)
        first(ks) = ielt
!       write (6,*) 'ielt,k,ks,next',ielt,k,ks,next(ielt)
      end do
! reset keep%varflag
      keep%varflag(1:n) = 0
    end if
! Note : in row entry case, we count only the generated elements

    mx_nchild = 0
    do 200 jdum = 1, n
      ivar = abs(perm(jdum))
!     write (6,*) 'jdum,ivar,order(ivar)',jdum,ivar,order(ivar)
      if (ivar == 0) cycle
! Check whether we have already dealt with this variable.
      if (order(ivar) == 0) cycle
      is = keep%svar(ivar)
!    write (6,'(a,9i4)') 'jdum,ivar,is,keep%nsup', jdum,ivar,is,keep%nsup

      ie = first_free
      if (ie < 0) then
! New node number needed.
        ie = keep%tnode + 1
        keep%tnode = ie
! Initialize node as empty
        keep%size(ie) = 0
        keep%tree(ie)%nelim = 0
        keep%ifile(ie) = 1
        keep%rfile(ie) = 1
      else
        first_free = next(ie)
      end if

! Check array child is large enough to hold children of ie
      j = -order(ivar)
      if (j < 0) then
        nchild = nels(ivar)
      else
! 2x2 pivot. Need to consider both entries in pivot.
        kvar = perm(j+1)
        nchild = nels(ivar) + nels(kvar)
      end if
      if (size(child) < nchild) then
        deallocate (child,stat=st)
        allocate (child(1:nchild),stat=st)
        if (st /= 0) return
      end if

! Construct the child list for ie.
      ielt = first(is)
      do i = 1, nels(ivar)
        child(i) = ielt
        ielt = next(ielt)
      end do
      ks = is
      if (j > 0) then
! 2x2 pivot.
        ks = keep%svar(kvar)
        ielt = first(ks)
        do i = 1, nels(kvar)
          child(nels(ivar)+i) = ielt
          ielt = next(ielt)
        end do
      end if
!     write (6,*) 'children:',child(1:nchild)

! nelim holds number of variables that can be eliminated at ie.
! The indices of these variables are placed in the array list from the front
! (keeping 2x2 pivots together).
! The indices of the other variables in the front are placed in list from
! the back; ifs holds the position of the first such variable.
      nelim = 0
      ifs = llist + 1

      if (.not.keep%element_input) then
! Row entry.
! Read the variable list for row ivar from the main integer
! superfile into varlist. Check diagonal present.
        l = 0
        nvar = abs(keep%size(ivar))
        loc = keep%ifile(ivar)
        call read_varlist(nvar)
        if (info%flag < 0) return
        do ii = 1, nvar
          jvar = varlist(ii)
          if (order(jvar) == 0) cycle
          if (jvar == ivar) l = 1
! jvar not yet eliminated.
          js = keep%svar(jvar)
          if (js == is .or. js == ks) then
! jvar belongs to the same supervariable as ivar or kvar and so
! can be eliminated.
            nelim = nelim + 1
            list(nelim) = jvar
            keep%varflag(jvar) = nelim
          else
! Put jvar at the end of list as not ready for elimination
            ifs = ifs - 1
            list(ifs) = jvar
            keep%varflag(jvar) = ifs
          end if
        end do
        if (l == 0) then
! add in diagonal
          nelim = nelim + 1
          list(nelim) = ivar
          keep%varflag(ivar) = nelim
        end if

        l = 0
        if (ks /= is) then
! Read in row kvar and then loop over entries in row kvar.
! Check diagonal present.
!         write (6,*) 'read row',kvar
          nvar = abs(keep%size(kvar))
          loc = keep%ifile(kvar)
          call read_varlist(nvar)
          if (info%flag < 0) return
          do ii = 1, nvar
            jvar = varlist(ii)
            if (order(jvar) == 0) cycle
! cycle if jvar already been encountered when looping over entries in row ivar.
            if (keep%varflag(jvar) > 0) cycle
            if (jvar == kvar) l = 1
            js = keep%svar(jvar)
            if (js == is .or. js == ks) then
              nelim = nelim + 1
              list(nelim) = jvar
              keep%varflag(jvar) = nelim
            else
              ifs = ifs - 1
              list(ifs) = jvar
              keep%varflag(jvar) = ifs
            end if
          end do
          if (l == 0) then
! Diagonal not present in row kvar. Must add if not
! already encountered kvar when dealing with row ivar.
            if (keep%varflag(kvar) <= 0) then
              nelim = nelim + 1
              list(nelim) = kvar
              keep%varflag(kvar) = nelim
            end if
          end if
        end if
      end if

! Now loop over the children of ie (which are all generated elements
! in the row entry case), reading in the variables
! for each child that are not yet eliminated and adding them into list.
        do i = 1, nchild
          ielt = child(i)
          nelim_i = 0
          if (ielt > keep%nelt) nelim_i = keep%tree(ielt)%nelim
          nvar = abs(keep%size(ielt)) - nelim_i
          loc = keep%ifile(ielt) + nelim_i
          call read_varlist(nvar)
          if (info%flag < 0) return
! varlist holds the variables in element ielt.
! Set negative flags so that we can quickly check if a partner is present
            do ii = 1,nvar
              jvar = varlist(ii)
              if (keep%varflag(jvar) == 0) then
                 keep%varflag(jvar) = -(llist+1)
              else
                 keep%varflag(jvar) = -keep%varflag(jvar)
              end if
            end do
! Merge varlist into list.
            do ii = 1, nvar
              jvar = varlist(ii)
              if (keep%varflag(jvar) == -(llist+1)) then
! First encounter of jvar for element ie.
                js = keep%svar(jvar)
                if (js == is) then
! jvar belongs to the same supervariable as ivar and so can be eliminated
                  nelim = nelim + 1
                  list(nelim) = jvar
                  keep%varflag(jvar) = -nelim
                else
! Put jvar at the end of list
                  ifs = ifs - 1
                  list(ifs) = jvar
                  keep%varflag(jvar) = -ifs
                end if
              end if
! Decrement nels(jvar)
              nels(jvar) = nels(jvar) - 1
!             write (6,*) 'jvar,nels(jvar)',jvar,nels(jvar)
              kk = -order(jvar)
              if (kk < 0) cycle
! jvar is part of 2x2 pivot. Check whether jvar is 1st or 2nd in the 2x2 pivot.
! If it is the 2nd entry and its partner jjvar is in ielt, do not
! decrement nels(jvar).
              if (perm(kk) < 0) cycle
              jjvar = -perm(kk-1)
              if (keep%varflag(jjvar) < 0) nels(jvar) = nels(jvar) + 1
            end do
! restore signs on keep%varflag
            do ii = 1,nvar
              jvar = varlist(ii)
              keep%varflag(jvar) = abs(keep%varflag(jvar))
            end do
        end do

! write (6,*) nelim,ifs,keep%varflag(2)
! write (6,*) list(1:nelim)
! write (6,*) list(ifs:llist)

      if (keep%element_input) then
! Element entry.
! Check whether any of the variables in list(ifs:llist) are fully
! summed (that is, nels(jvar) = 0). If they are, move to front of array list.
        k1 = ifs
        do k = k1,llist
          jvar = list(k)
          if (nels(jvar) == 0) then
! jvar can be eliminated so move into front portion of list
            nelim = nelim + 1
            list(nelim) = jvar
            keep%varflag(jvar) = nelim
! make sure we leave no gaps in list
            if (k /= ifs .and. ifs < llist) then
              kvar = list(ifs)
              list(k) = kvar
              keep%varflag(kvar) = k
            end if
            ifs = ifs + 1
          end if
        end do

      end if

!         write (6,*) 'nels',nels(1:n)
!         write (6,*) 'nelim,list 1:',nelim,list(1:nelim)
!         write (6,*) 'ifs,list 2:',ifs,list(ifs:llist)
!         write (6,*) 'perm:',perm(1:n)

! Check that if a variable that is part of a 2x2 pivot is marked for
! elimination then its partner is also ready for elimination
        ii = 0
        do k = 1,nelim
          jvar = list(k)
          j = -order(jvar)
          if (j < 0) cycle
          if (perm(j) > 0) then
            kvar = -perm(j-1) ! second entry in 2x2
          else
            kvar = perm(j+1)
          end if
          if (keep%varflag(kvar) == 0 .or. keep%varflag(kvar) >= ifs) then
            ii = ii + 1
            list(k) = -jvar
          end if
        end do

        if (ii > 0) then
! Move those that could not be eliminated to the end section of list
          l = 0
          do k = 1,nelim
            jvar = list(k)
!          write (6,*) 'k,jvar,l,ifs',k,jvar,l,ifs
            if (k == ifs) exit
            if (jvar == 0) exit
            list(k) = 0
            if (jvar > 0) then
              l = l + 1
              list(l) = jvar
              keep%varflag(jvar) = l
            else
              ifs = ifs - 1
              if (ifs > nelim) then
                list(ifs) = -jvar
                keep%varflag(-jvar) = ifs
              else
                k1 = list(ifs)
                do while (k1 < 0)
                  list(ifs) = -jvar
                  keep%varflag(-jvar) = ifs
                  ifs = ifs - 1
                  k1 = list(ifs)
                end do
                list(ifs) = -jvar
                keep%varflag(-jvar) = ifs
                if (k1 == 0) exit
                l = l + 1
                list(l) = k1
                keep%varflag(k1) = l
              end if
            end if
          end do
          nelim = l
        end if

!       write (6,*) 'nelim,list 1:',nelim,list(1:nelim)
!       write (6,*) 'order',order(1:n)
! We need to reunite any pairs that have become separated within list(1:nelim).
! We flag the first entry in the pair with a negative sign.
         l = 1
         do i = 1,nelim
           if (l > nelim) exit
           jvar = list(l)
           j = -order(jvar)
           if (j < 0) then
             l = l + 1; cycle
           end if
! jvar is part of a pair. Let its partner be kvar.
           k1 = list(l+1)
           if (perm(j) < 0) then
! jvar is the first in the pair.
             kvar = perm(j+1)
             k = keep%varflag(kvar)
             list(l) = -jvar
             list(l+1) = kvar
             list(k) = k1
             keep%varflag(jvar) = l
             keep%varflag(kvar) = l+1
             keep%varflag(k1) = k
           else
! jvar is the second in the pair.
             kvar = -perm(j-1)
             k = keep%varflag(kvar)
             list(l) = -kvar
             list(l+1) = jvar
             if (k > l+1) list(k) = k1
             keep%varflag(kvar) = l
             keep%varflag(jvar) = l+1
             if (k > l+1) keep%varflag(k1) = k
           end if
           l = l + 2
           info%ntwo = info%ntwo + 1
         end do

  !      write (6,*) 'At node= ', ie, 'nelim= ',nelim, &
  !    ' list: ',list(1:nelim),list(ifs:llist)
! Set keep%varflag to -ie for variables that are eliminated
! and set > 0 for those that are not eliminated
      do i = 1,nelim
         jvar = abs(list(i))
         order(jvar) = 0
         keep%varflag(jvar) = -ie
      end do
      do i = ifs,llist
         jvar = list(i)
         keep%varflag(jvar) = i
      end do
! nfs is the number of non fully summed variables
      nfs = llist - ifs + 1

      if (nfs > 0) then
! Link new element into list. Loop over the uneliminated variables belonging
! to ie to find the one that occurs earliest in the pivot sequence.
        jj = n + 1
        do i = ifs,llist
          jvar = list(i)
          nels(jvar) = nels(jvar) + 1
          kk = -order(jvar)
          belong = .false.
          if (kk > 0) then
! jvar is part of 2x2 pivot. Check whether it is 1st or 2nd in the 2x2 pivot.
! If it is the 2nd entry, check to see if its partner kvar is in list
            if (perm(kk) > 0) then
              kvar = -perm(kk-1)
              if (keep%varflag(kvar) > 0) then
                nels(jvar) = nels(jvar) - 1
                belong = .true.
              end if
            end if
          end if
          if (abs(order(jvar)) < jj) then
            k = jvar
            if (belong) k = kvar
            jj = abs(order(k))
          end if
        end do

        ks = keep%svar(k)
        next(ie) = first(ks)
        first(ks) = ie

      else
! Node is a root. Link with other roots
        next(ie) = first_root
        first_root = ie
        nroot = nroot + 1
      end if

! reset keep%varflag to zero for variables that are not eliminated
      do i = ifs,llist
         jvar = list(i)
         keep%varflag(jvar) = 0
      end do

! Merge parent ie and child nodes if the list of uneliminated variables
! at the child is the same as the list of variables at the parent
! or if both involve fewer than nemin eliminations.
      k = nchild
      do 80 i = 1, k
        ielt = child(i)
! Cycle if ielt is a leaf node
        if (ielt <= keep%nelt) cycle
! Skip to next child if merging is not needed.
        nelim_i = keep%tree(ielt)%nelim

        if (nelim+nfs > keep%size(ielt)-nelim_i .and. &
           (nelim >= nemin .or. nelim_i >= nemin)) cycle

! Merge. Read in variable list. Only need the first nelim_i entries for ielt.
! Read directly into list
        if (nelim_i /= 0) then
          loc = keep%ifile(ielt)
          call MA77_read_integer(keep%index(1),keep%imain,loc,nelim_i,    &
               list(nelim+1:nelim+nelim_i),flag,keep%idata,-1)
          if (loc+keep%size(ielt) == keep%ifree) &
! store the largest integer storage needed so far
            keep%ifree = loc
          if (flag < 0) then
            info%flag = flag; return
          end if
          nelim = nelim + nelim_i
        end if

! The children of ielt become children of ie.
        nc = size(keep%tree(ielt)%child)
! ie may have more children than before so check size
! of array child. If necessary, allocate larger array child.
        if (nchild+nc-1 > size(child)) then
          call reallocate_child(nchild+nc-1)
          if (info%flag < 0) return
        end if

! Loop over children of ielt, making them children of ie
        if (nc > 0) then
! The first child replaces ielt in the list of children for ie
          l = i
          do j = 1, nc
            jelt = keep%tree(ielt)%child(j)
            child(l) = jelt
            l = nchild + j
          end do
          nchild = nchild + nc - 1
        else
! Flag that the child has gone.
          child(i) = 0
        end if

! Remove node ielt from keep%tree by deallocating its list of children
! and adding it to the list of free nodes.
        deallocate (keep%tree(ielt)%child,stat=st)
        next(ielt) = first_free
        first_free = ielt
        keep%size(ielt) = 0
        keep%tree(ielt)%nelim = 0
        keep%ifile(ielt) = 1
        keep%rfile(ielt) = 1

! Finished with child ielt
80    continue

! Flush out dummy children
      k = nchild
      nchild = 0
      do i = 1, k
        if (child(i) /= 0) then
          nchild = nchild + 1
          child(nchild) = child(i)
        end if
      end do
      mx_nchild = max(mx_nchild,nchild)

! Copy list of children into tree(ie)%child
      allocate (keep%tree(ie)%child(1:nchild),stat=st)
    ! write (6,*) 'node',ie,' has children',child(1:nchild)
      keep%tree(ie)%child(1:nchild) = child(1:nchild)

! Move list of uneliminated variables forward
      k = nelim
      do i = ifs,llist
         k = k + 1
         list(k) = list(i)
      end do

      keep%size(ie) = k
      keep%maxlen = max(keep%maxlen,k)
      keep%maxfront = max(keep%maxfront,k)
!     write (*,*) 'ie,k,maxlen',ie,k,keep%maxlen
      keep%tree(ie)%nelim = nelim
      keep%maxelim = max(nelim,keep%maxelim)

! Store variable list for ie in main integer superfile.
      loc = keep%ifree
!       write (6,*) 'write',k,' integers',' from position',loc
!       write (6,*) 'ie,nelim,list',ie,nelim,list(1:k)

      call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
           loc,k,list,flag,keep%idata,-1,keep%maxstore,keep%used)
      if (flag < 0) then
        info%flag = flag
        return
      end if

! Store position in the integer superfile of this variable list
      keep%ifile(ie) = loc
! Move pointer for next free location in the integer superfile
      keep%ifree = loc + k
      keep%mx_ifree = max(keep%mx_ifree,keep%ifree)

!!! end of main tree loop
200 continue
    info%tree_nodes = keep%tnode - keep%nelt


   end subroutine construct_tree2

!**********************************

!!! End of analyse

  end subroutine MA77_analyse_double

!****************************************************************************

  subroutine MA77_input_reals_double(index,length,reals,keep,control,info)

! This subroutine is called by the user to specify the reals
! for an element/row (reverse communication interface).
! For element entry, the lower triangular part of the element matrix
! must be input by columns in packed form.
! If there are no duplicates/out-of-range entries in the
! element/row and if the sequence of nonzeros
! is very large, it may be broken into parts and input by
! consecutive calls to MA77_input_reals.
! For row entry, the user must input the nonzeros
! in the row (both upper and lower triangular entries are needed).
! Again, a number of consecutive calls may be used for one row.

! The integer data for an element/row must be entered
! before the reals for the elt/row (but all the calls to MA77_input_vars
! do not have to be completed before the calls to MA77_input_reals
! are started).
! If values are entered for an element/row that has already
! been entered, the new values overwrite the old (this allows
! for subsequent factorizations of further matrices
! with same pattern but different numerical values).

    integer(short), intent (in) :: index ! must hold the index of
!            incoming element/row. If out of range, the call is ignored.
    integer(short), intent (in) :: length ! must hold the number of nonzeros
!            being input. Must be  >=  0.
    real (wp), intent (in) :: reals(length) ! must contain the nonzeros
!            being input
! For details of keep, control, info : see derived type description
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info

! Local scalars
    integer(short) :: flag ! local error flag
    integer(short) :: i,ii
    integer(short) :: inelrs
    integer(short) :: j,jj
    integer(short) :: k,kk
    integer(short) :: len
    integer(short) :: nout ! unit for printing errors
    integer(short) :: nvar ! no. of variables in incoming element/row
!              (compressed list)
    integer(short) :: nvar_user ! no. of variables in user-supplied list for
!              for incoming element/row (before duplicated/out-of-range
!              entries removed)
    integer(short) :: st
    integer(long) :: loc,locw

! Possible error returns:
! -1   Allocation error (only possible if duplicates and/or out-of-range
!      entries)
! -3   Error in sequence of calls to routines in package
! -5   Error from of01_write (error in Fortran inquire)
! -6   Error from of01_write (error in Fortran read)
! -7   Error from of01_write (error in Fortran open)
! -10  MA77_input_vars was not called
!      for the element/row index. No action taken.
! -11  Problem found to be singular or unexpectedly not pos.def.
! -14  Data for previous element/row incomplete. No action taken.
! -15  Error from of01_write (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!      for all path)
! -17  out-of-range and/or duplicated indices entered but
!      user has not entered all the reals for the current elt/row
!      in a single call to input_reals. No action taken.
! -19  length < 0. No action taken.
! -32  Too many reals input for current element/row. No action taken.

    if (keep%n == 0) return

! Perform appropriate printing
    if (control%print_level >= 2 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a,2(/a,i12))') &
     ' Entering MA77_input_reals with:', &
     ' index               Index of element/row                = ',index, &
     ' length              Number of reals being input         = ',length
      i = min(length,10)
      write (control%unit_diagnostics,'(a,2(/5es12.4))') &
     ' Input reals: ',reals(1:i)
      if (i < length) write (control%unit_diagnostics,'(a)') '  . . . . . .'
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_input_reals'

! Remember that calls to this routine may follow an earlier factorization
! and so we do not want to retain a flag for singularity set on the earlier
! call. 
    if(keep%flag == 4) keep%flag = 0

    flag = 0; info%flag = keep%flag

! Check status parameter (remember that input of integer and real
! data can be interleaved and user can choose not to input reals
! until after MA77_analyse has been called).
! Remember also that call could follow a factorization.
    if (keep%status  < 1) then
      info%flag = -3
! Check input data
    else if (index < 1 .or. index > keep%nelt) then
      return
    else if (keep%ifile(index) == 0) then
! no integer data has been entered for an element/row of this index.
      info%flag = -10
    else if (length  < 0) then
      info%flag = -19
    end if

    if (info%flag < 0) then
      keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      return
    end if

! Immediate return if keep%size(index) = 0 (this can happen only
! if nvar = 0 on call to MA77_input_vars)
    if (keep%size(index) == 0) return

! We will need some arrays to squeeze out duplicates etc
! Note: input of variables and reals can be interleaved so at this stage
! keep%mvar holds largest value of nvar so far.
    if (keep%size(index) < 0) then
      st = 0
      if (keep%element_input) then
        len = (keep%mvar*keep%mvar+keep%mvar)/2
        if (.not. allocated(keep%aelt)) then
           allocate (keep%aelt(len),stat=st)
        else if (size(keep%aelt) < len) then
           deallocate (keep%aelt,stat=st)
           allocate (keep%aelt(len),stat=st)
        end if
      else
        len = keep%mvar
        if (.not. allocated(keep%arow)) then
           allocate (keep%arow(len),stat=st)
        else if (size(keep%arow) < len) then
           deallocate (keep%arow,stat=st)
           allocate (keep%arow(len),stat=st)
        end if
      end if
      if (st /= 0) then
        info%stat = st
        info%flag = -1; keep%status = -1
        if (nout >= 0) &
          call MA77_print_iflag(keep,nout,info%flag,ie=index,st=info%stat)
        return
      end if
    end if

! keep%inelrs is the number of reals still needed for the current
! element/row. If it is equal to zero on input
! then no reals have yet been input for the current element/row.
! if keep%inelrs  >  0,  keep%inelrn holds the index of the
! element/row that was last input.
!   write (6,*) 'index,info%flag,inelrs,length,keep%inelrn',&
!   index,info%flag,keep%inelrs,length,keep%inelrn
    inelrs = keep%inelrs
    if (inelrs > 0) then
      if (index /= keep%inelrn) then
! data for previous element/row is not complete
        info%flag = -14
        call MA77_print_iflag(keep,nout,info%flag,ie=index)
        return
      end if
    else if (index == keep%inelrn .and. length > 0) then
! Too many reals have been supplied for element/row
      info%flag = -32
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      return
    else
! no reals have yet been input for this element/row.
! Set inelrs to the expected no. of entries
      if (keep%size(index) > 0) then
        keep%inelrn = index
        nvar = keep%size(index)
        inelrs = nvar
        if (keep%element_input) inelrs = (nvar*nvar + nvar)/2
      else
! some duplicates/out-of-range variables were entered for this element/row.
! Read compressed list and length of original user-supplied list (nvar_user)
        nvar = -keep%size(index)
        loc = keep%ifile(index)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar+1,  &
             keep%clist,flag,keep%idata,-1)
        if (flag < 0) go to 50
        nvar_user = keep%clist(nvar+1)
!       write (6,*) 'clist,nvar_user',keep%clist(1:nvar),nvar_user
        inelrs = nvar_user
        if (keep%element_input) inelrs = (nvar_user*nvar_user + nvar_user)/2
! We require the user to enter the whole element/row using a single call
! If supplied incorrect number, the user can correct length and recall
        if (inelrs > length) then
          info%flag = -17
          call MA77_print_iflag(keep,nout,info%flag,ie=index)
          return
        end if
      end if
    end if

! If length = 0, nothing to be done
    if (length == 0) return

! check whether the user has supplied too many reals.
! If too many supplied, the user can correct length and recall
! (keep%inelrs has not been altered)
    if (length > inelrs) then
      info%flag = -32
      call MA77_print_iflag(keep,nout,info%flag,ie=index)
      return
    end if

! Store the reals in main real file.
! If we already have real data for an earlier version
! of the current element/row, we overwrite it with latest real data
    if (keep%rfile(index) > 0 .and. keep%inelrs == 0) then
      locw = keep%rfile(index) - inelrs
    else
      locw = keep%posfac
! move pointer for next free location in real file
      if (keep%size(index) > 0) then
         keep%posfac = keep%posfac + length
      else
         len = nvar
         if (keep%element_input) len = (nvar*nvar + nvar)/2
         keep%posfac = keep%posfac + len
      end if
    end if

! If no duplicates/out-of-range entries in this element/row, we can write
! directly to the main real file
    if (keep%size(index) > 0) then
      call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain,&
           locw,length,reals,flag,keep%rdata,-1,keep%maxstore,keep%used)
    else
! We have to squeeze out duplicates/out-of-range entries.
! Start by reading in the mapping that was stored by MA77_input_vars.
      loc = keep%ifile(index) + nvar + 1
      call MA77_read_integer(keep%index(1),keep%imain,loc,nvar_user,&
           keep%map(1:nvar_user),flag,keep%idata,-1)
!     write (6,*) 'loc,map',loc,keep%map(1:nvar_user)
      if (flag < 0) go to 50
      if (keep%element_input) then
! Set ip(jj) to point to start of col jj in array aelt.
        keep%iptr(1) = 1
        do jj = 2,nvar
          keep%iptr(jj) = keep%iptr(jj-1) + nvar - jj + 2
        end do
        k = 1
        keep%aelt(1:len) = zero
! Loop over the columns of the incoming element
        do j = 1,nvar_user
          jj = keep%map(j)
! Skip if column corresponds to an out-of-range index
          if (jj == 0) then
            k = k + (nvar_user - j + 1)
            cycle
          end if
! Loop over the rows in column j.
          do i = j,nvar_user
            ii = keep%map(i)
            if ( ii /= 0) then
              if (ii >= jj) then
                kk = keep%iptr(jj) + ii - jj
              else
                kk = keep%iptr(ii) + jj - ii
              end if
              keep%aelt(kk) = keep%aelt(kk) + reals(k)
            end if
            k = k + 1
          end do
        end do
        call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain,&
             locw,len,keep%aelt,flag,keep%rdata,-1,keep%maxstore,keep%used)
      else
        keep%arow(1:len) = zero
        do i = 1,nvar_user
          j = keep%map(i)
          if (j == 0) cycle
          keep%arow(j) = keep%arow(j) + reals(i)
        end do
        call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain,&
             locw,len,keep%arow,flag,keep%rdata,-1,keep%maxstore,keep%used)
      end if
    end if
    if (flag < 0) go to 60

! Set keep%inelrs to hold number of reals still required for current
! element/row
    keep%inelrs = inelrs - length

! Set rfile(index) to point to the position in the real file
! after the current element/row (it points to the end of the
! elt/row list rather than to the start to allow us overwrite the data
! with new real values in the case when the elt/row is input
! using more than one call ... if more than one call is used
! and we want to overwrite existing data, we cannot do this
! if we only know the position of the start of the list for the elt/equ)
    keep%rfile(index) = locw + length
    if (keep%size(index) < 0) keep%rfile(index) = locw + len

! keep%posfac holds first free location in main real  
! superfile where we can write factor

    if (control%print_level >= 2 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a,(/a,i12))') &
     ' Leaving MA77_input_reals with:', &
     ' Error flag                      info%flag             = ', info%flag
    end if
    return

 50 info%iostat = keep%idata%iostat

 60 info%iostat = keep%rdata%iostat
    info%flag = flag; keep%status = -1
    call MA77_print_iflag(keep,nout,info%flag,ie=index,ios=info%iostat)

  end subroutine MA77_input_reals_double

!****************************************************************************

  subroutine MA77_factor_double(pos_def,keep,control,info,scale)

! Factorisation phase.

    logical :: pos_def ! Must be set to .true. if the problem is known to
!            positive definite and to .false. otherwise.
!            if pos_def = .true., no pivoting is performed
! For details of keep, control, info : see derived type description
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info

    real(wp), intent(in), optional :: scale(:)  ! if present, must hold
!     row and column scaling factors

    real(wp) :: x(0,0)

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -11   Problem found to be singular or unexpectedly not pos.def.
! -14   data for last element/row passed to MA77_input_reals was incomplete.
! -15   Error from of01 (error in Fortran write)
! -22   For one or more elements/rows, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.
! -29   pos_def = .true. but 2x2 pivots supplied to analyze
! -30   front size is too large
! -38   error in static
! -39   size(scale) too small
! -40 IEEE infinities found in the reduced matrix

!   write(*,*) 'entering MA77_factor with keep%status=', keep%status

    keep%name = 'MA77_factor'

    if (present(scale)) then
      call MA77_factor_solve_double(pos_def,keep,control,info,0,0,x,scale)
    else
       call MA77_factor_solve_double(pos_def,keep,control,info,0,0,x)
    end if

   end subroutine MA77_factor_double

!****************************************************************************

  subroutine MA77_factor_solve_double(pos_def,keep,control,info,nrhs, &
     lx,x,scale)

! Factorisation phase. If lx .ge. n and nrhs > 0, also solve.

    integer(short), parameter :: nb_default = 150
    logical, intent (in) :: pos_def ! Must be set to .true. if
!            the problem is known to
!            positive definite and to .false. otherwise.
!            if pos_def = .true., no pivoting is performed
! For details of keep, control, info : see derived type description
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info

    real(wp), intent(in), optional :: scale(:)  ! if present, must hold
!     row and column scaling factors

    integer(short) :: lx ! First extent of x. must be at least n.
    integer(short) :: nrhs ! number right-hand sides. must be at least 1
!          if solve required
    real (wp), intent(inout) :: x(lx,nrhs) ! if nrhs > 0
!          must hold right-hand sides on entry. On exit, holds solution.

! The following arrays are allocated by MA77_factor_solve and are
! then used as workspace by subroutine MA77_factorize
    real (wp), allocatable :: fa(:) ! frontal matrix
    real (wp), allocatable :: w54(:) ! array used by ma54 (pos. def. kernel)
    real (wp), allocatable :: w64(:) ! array used by ma64 (indefinite kernel)
    real (wp), allocatable :: reals(:) ! used to hold a single col.
!            of the frontal matrix. Also used for reading in
!            a user-supplied elt/equ.
    real (wp), allocatable :: xlocal(:,:) ! temporary dense vector for
!              performing forward subs (only needed if x present).
    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree
    integer(short), allocatable :: rlist(:), varlist(:) ! used for holding
!            the variables at nodes in the tree.
!            They are allocated to be of length keep%maxlen. This is
!            sufficient if there are no delayed cols. If insufficient,
!            reallocated within MA77_factorize.
    integer(short), allocatable :: map(:) ! used to map between variables in
!            a child and those in its parent.
    integer(short), allocatable :: pos(:) ! Used for mapping into the front
!            (only array of length n)
    integer(long), allocatable :: ip(:) ! Used to hold
!            positions of columns in the front.

    type(ma64_control) :: cntl64 ! control derived type for hsl_ma64
    type(ma64_info) :: info64 ! information derived type for hsl_ma64

    real (wp) :: multiplier ! set to max(1.0, control%multipler) (indef. case)
    logical :: lsolve ! set to true if solves are to follow the factorization

  integer(short) :: cdelay ! number of delayed pivots for children of code
  integer(short) :: cnvar ! number of (uneliminated) variables at child node
  integer(short) :: depth ! depth in tree
  integer(short) :: d1max ! max value of d(1,:) over children of node
  integer(short) :: dreals !
  integer(short) :: flag
  integer(short) :: i ! temporary variable
  integer(long)  :: ilong ! temporary variable
  integer(short) :: ic ! indicates which child of pnode node is
  integer(short) :: ii ! temporary variable
  integer(long)  :: inactive_fac ! points to first entry in factor storage
!   (= keep%posfac + n in pos def case; keep%posfac + 2n in indefinite case)
  integer(short) :: inelrs ! number of reals to be read in
  integer(long)  :: ipp ! points to start of a column in frontal matrix
  integer(long)  :: ip1 ! points to position of start of frontal mx for node
  integer(long)  :: ipc ! points to position of start of frontal
!       mx for child
  integer(long)  :: ipj2 ! temporary variable (used assemble54/assemble64)
  integer(long)  :: ipk  ! temporary variable (used assemble54/assemble64)
  integer(short) :: ir ! do loop variable
  integer(short) :: ivar ! a variable
  integer(short) :: j  ! temporary variable
  integer(long)  :: jlong  ! temporary variable
  integer(short) :: j1 ! temporary variable
  integer(short) :: j2 ! temporary variable
  integer(short) :: jc ! temporary variable
  integer(short) :: jj ! temporary variable
  integer(short) :: jcnode ! jc-th child node
  integer(short) :: jrow ! a row in the front
  integer(short) :: jvar ! a variable
  integer(short) :: k  ! temporary variable
  integer(short) :: k1 ! temporary variable (used in scaling)
  integer(short) :: k2 ! temporary variable (used in scaling)
  integer(short) :: kk ! temporary variable
  integer(long)  :: klong ! temporary variable
  integer(short) :: l  ! temporary variable
  integer(short) :: ld !
  integer(short) :: ldc ! number of delayed columns for child
  integer(long)  :: lfa ! size of array fa
  integer(long)  :: length ! total length of list to be read
  integer(long)  :: lnj,lnj1 ! used to compute flop count and factor entries
  integer(short) :: lreals ! number of reals to be written
!     Also length of array reals
  integer(long)  :: llreals !
  integer(long)  :: ln ! set to (nfront*(nfront+1))/2 during factorize
  integer(long)  :: llong
  integer(long)  :: loc1
  integer(long)  :: loc  ! location in real or integer superfile when reading
  integer(long)  :: locw ! location in real or integer superfile when writing
  integer(long)  :: loci ! location in integer superfile
  integer(long)  :: locr ! location in super-real file
  integer(long)  :: lq ! set by call to ma64_factor
  integer(short) :: lw54 ! length of array w54
  integer(short) :: lw64 ! length of array w64
  integer(short) :: lw
  integer(short) :: maxstore
  integer(short) :: maxlen ! In pos. def. case, holds copy of keep%maxlen
!           In indef. case, set to keep%maxlen*control%multiplier
  integer(short) :: mvar ! copy of keep%mvar
  integer(short) :: nb ! set to control%nb54 or control%nb64
!        (or to nb54_default or nb64_default if out of range)
  integer(short) :: nbi ! set to control%nbi
!        (or to nbi_default if out of range)
  integer(short) :: nfront ! order of frontal matrix to be factorized
  integer(long)  :: nfront_long ! order of frontal matrix to be factorized
  integer(short) :: node ! node in tree
  integer(short) :: nchild ! number of children
  integer(short) :: n ! set to keep%n
  integer(short) :: ncand ! number of variables that are candidates
!         for elimination at node
  integer(short) :: nelim_jc ! no. variables to be eliminated at child jcnode
  integer(short) :: nelim  ! number of eliminations performed at node
  integer(short) :: non_leaf ! number of children of node that are not
!         leaf nodes
  integer(short) :: nout ! output unit for errors
  integer(short) :: nout1 ! output unit for warnings
  integer(short) :: nroot ! number of roots (components) of the tree
  integer(short) :: nvar ! no. variables associated with node during analyse
  integer(long)  :: nvar_long ! no. variables associated with node 
!         during analyse
  integer(long)  :: posdiag ! points to end of storage for diagonal entries
!                    (initialised to keep%posfac)
  integer(short) :: pnode ! parent node in tree
  integer(short) :: root ! root node
  integer(long)  :: rtop ! points to the top of the main real stack
  integer(long)  :: rtopd ! points to the top of the real stack for
!           delayed pivots (only needed in indefinite case)
  integer(long)  :: size_fa ! size of frontal array fa
  integer(short) :: splitp ! split point at which space allocated for
!         generated frontal matrix
  integer(short) :: st ! stat parameter
  integer(long)  :: storage_indef
  integer(long)  :: temp(8) ! used to get counts of of01 reads/writes when
!           solve follows factorization

  integer(long), allocatable :: dpos(:) ! used to pass info.
!       on delayed pivots
  integer(short), allocatable :: delay(:,:) ! Set so that
! delay(1,node) holds nfront-nelim where nfront is frontsize
!      and nelim is no. of eliminations performed
! delay(2,node) holds number of delayed pivots ie ncand-nelim
!      where ncand is no. of candidate pivots (no. of fully summed variables)

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -11   Problem found to be singular or unexpectedly not pos.def.
! -14   data for last element/row passed to MA77_input_reals was incomplete.
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -22   For one or more elements/rows, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.
! -24   lx or nrhs out of range
! -29   pos_def = .true. but 2x2 pivots supplied to analyze
! -30   front size is too large
! -38   error in static
! -39   size(scale) too small
! -40 IEEE infinities found in the reduced matrix

!  +4   Problem found to be singular (control%action = .true.)

! Temporary variables for timing
!    real :: t1, t2, tfactor, t54, st1, st2, stfactor, st54, tmerge, stmerge
!    tfactor = 0.0; t54 = 0.0; stfactor = 0.0; st54 = 0.0
!    tmerge = 0.0; stmerge = 0.0

    lsolve = .true.
    if (lx == 0 .and. nrhs == 0) lsolve = .false.

!   write(*,*) 'entering MA77_factor_solve with keep%status=', keep%status
! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      if (pos_def) then
       if (lsolve) write (control%unit_diagnostics,'(//a)') &
       ' Entering MA77_factor_solve with pos_def = .true. and :'
       if (.not. lsolve) write (control%unit_diagnostics,'(//a)') &
       ' Entering MA77_factor with pos_def = .true. and :'
        write (control%unit_diagnostics,'(a,5(/a,i12))') &
       ' control parameters (control%) :', &
       ' print_level         Level of diagnostic printing           = ', &
         control%print_level,  &
       ' unit_diagnostics    Unit for diagnostics                   = ', &
         control%unit_diagnostics, &
       ' unit_error          Unit for errors                        = ', &
         control%unit_error,   &
       ' unit_warning        Unit for warnings                      = ', &
         control%unit_warning, &
       ' nb54                Block size for HSL_MA54                = ', &
         control%nb54
      else
       write (control%unit_diagnostics,'(//a)') &
       ' Entering MA77_factor with pos_def = .false. and :'
        write (control%unit_diagnostics,&
         '(a,6(/a,i12),6(/a,es12.4))') &
       ' control parameters (control%) :', &
       ' print_level         Level of diagnostic printing           = ', &
         control%print_level,      &
       ' unit_diagnostics    Unit for diagnostics                   = ', &
         control%unit_diagnostics, &
       ' unit_error          Unit for errors                        = ', &
         control%unit_error,       &
       ' unit_warning        Unit for warnings                      = ', &
         control%unit_warning,     &
       ' nb64                Block size for HSL_MA64                = ', &
         control%nb64,             &
       ' nbi                 Inner block size for HSL_MA64          = ', &
         control%nbi,              &
       ' small               Small pivot size                       = ', &
         control%small,            &
       ' static              Static pivoting control                = ', &
         control%static,           &
       ' u                   Initial relative pivot tolerance       = ', &
         control%u,                &
       ' umin                Minimum relative pivot tolerance       = ', &
         control%umin,             &
       ' multiplier          Multiplier for increasing array sizes  = ', &
         control%multiplier,       &
       ' storage_indef       Initial size for in-core array         = ', &
         real(control%storage_indef)
      end if
      if (lsolve) then
        write (control%unit_diagnostics,'(a,i12/a,i12)') &
       ' lx                                                         = ', &
         lx, &
       ' nrhs                                                       = ', &
         nrhs
      end if
    end if

    keep%scale = 0
    if (present(scale)) then
      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
      write (control%unit_diagnostics,'(a)') &
     ' Scaling factors supplied. '
      keep%scale = 1
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    nout1 = control%unit_warning
    if (control%print_level < 0) nout1 = -1
    if (lsolve) keep%name = 'MA77_factor_solve'

    info%flag = 0

    info%nio_read(1) = keep%idata%nio_read
    info%nio_read(2) = keep%rdata%nio_read
    info%nio_write(1) = keep%idata%nio_write
    info%nio_write(2) = keep%rdata%nio_write
    info%nwd_read(1) = keep%idata%nwd_read
    info%nwd_read(2) = keep%rdata%nwd_read
    info%nwd_write(1) = keep%idata%nwd_write
    info%nwd_write(2) = keep%rdata%nwd_write

! Check status parameter (status = 3 is possible if restart has
! been called)
    if (keep%status < 2) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (keep%n == 0) then
      keep%status = 3; return
    end if

! check the last call to MA77_input_reals was complete.
    if (keep%inelrs > 0) then
! data for last element/row not complete
      info%flag = -14
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    keep%pos_def = pos_def
    if (keep%pos_def) then
      nb = control%nb54
! check nb is valid ... if not, we will use our recommended value
      if (nb < 1) nb = nb54_default
! take copy for use in solve phase
      keep%nb = nb
      if (keep%ntwo > 0) then
! 2x2 pivots supplied
        info%flag = -29
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if

    else

! Set up controls for HSL_MA64
      nb = control%nb64
      nbi = control%nbi
! check nbi is valid ... if not, we will use our recommended value
      if (nbi <= 1) nbi = nbi_default

! check nb is valid ... if not, we will use our recommended value
      if (nb < 1) nb = nb64_default
! nb has to be a multiple of nbi. If necessary, reset to that it is.
      if (nb/nbi*nbi /= nb) nb = max(nbi,nb/nbi*nbi)
! take copy for use in solve phase
      keep%nb = nb
      keep%nbi = nbi
      cntl64%small = control%small

      cntl64%twos = .false.
      if (keep%ntwo > 0) cntl64%twos = .true.

! Check control%static
      if (control%static < control%small .and. control%static /= zero) then
        info%flag = -38
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if
      cntl64%static = control%static

! check threshold ... if out-of-range, use the default
      cntl64%u = control%u
      if (cntl64%u < zero .or. cntl64%u > one) cntl64%u = 0.01

      cntl64%umin = control%umin
      if (cntl64%umin < zero) cntl64%umin = zero
      if (cntl64%umin > cntl64%u) cntl64%umin = cntl64%u

! Check to see if the user is working in-core (note: if we are
! restarting the computation, keep%rwdelay is already allocated)
      if (.not. allocated(keep%rwdelay)) then
        maxstore = keep%maxstore
        if (maxstore > 0_long) then
          storage_indef = control%storage_indef
! user working in-core
          if (storage_indef == 0_long) then
            llong = max(1_long,maxstore/200_long)
            storage_indef = llong*4
          end if
          if (storage_indef > 0_long .and. &
            keep%used + 2*storage_indef <= maxstore) then
            deallocate (keep%rwdelay,stat=st)
            allocate (keep%rwdelay(storage_indef),stat=st)
            if (st == 0) then
              keep%index(4) = -keep%index(4)
              keep%used = keep%used + 2*storage_indef
              keep%size_rwdelay = storage_indef
            end if
          end if
        end if
      end if
    end if

    n = keep%n

    if (lsolve) then
      if (lx < n .or. nrhs < 1) then
        info%flag = -24
        keep%status = -1
        call MA77_print_iflag(keep,nout,info%flag)
        if (nout >= 0) then
          if (lx < n) write (nout,'(a,i8,a,i8)') &
      ' Increase lx from ', lx, ' to at least ', n
          if (nrhs < 1) write (nout,'(a,i8,a,i8)') &
      ' nrhs must be at least 1. nrhs = ', nrhs
        end if
        return
      end if
    end if

    if (present(scale)) then
      if (size(scale) < n) then
        info%flag = -39
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if
      if (lsolve) then
! scale right-hand sides
        do i = 1,n
          x(i,1:nrhs) = scale(i)*x(i,1:nrhs)
        end do
      end if
    end if


! Deallocate arrays used by analyse but no longer needed.
    deallocate (keep%clist,stat=st)
    deallocate (keep%arow,stat=st)
    deallocate (keep%aelt,stat=st)
    deallocate (keep%iptr,stat=st)

! Allocate arrays needed by factorize
    multiplier = one
    if (.not. keep%pos_def) multiplier = max(control%multiplier, one)
    mvar = keep%mvar

!    write (6,*) 'keep%index ',keep%index(1:4)
    deallocate (fa,stat=st)
    nfront_long = keep%maxfront
    if (keep%pos_def) then
      lfa = (nfront_long*(nfront_long+1))/2
      if (lfa > keep%lup) info%flag = -30
      maxlen = keep%maxlen
    else
      nfront_long = nfront_long*multiplier
      lfa = (nfront_long*(nfront_long+nb+1)/2)
      if (lfa > keep%lup) info%flag = -30
      maxlen = keep%maxlen*multiplier
    end if
    if (info%flag == -30) then
      call MA77_print_iflag(keep,nout,info%flag)
      go to 25
    end if

    allocate (fa(lfa),stat=st)

    if (st /= 0) then
! if arrays are being used, we try and switch to using files
      if (all(keep%index(1:4) > 0)) then
! files are being used so we have failed to allocate front
       info%flag = -30;  info%stat = st
       call MA77_print_iflag(keep,nout,info%flag)
       go to 25   
      else
! at least one array is being used.
! Try and write the integers to file
        if (keep%index(1) < 0) then
          keep%index(1) = -keep%index(1)
          loc1 = 1
          lw = int(keep%ifree,short) - 1
          call of01_write(keep%index(1),loc1,lw,keep%imain, &
             flag,keep%idata,lp=-1)
          if (flag < 0) then
            info%iostat = keep%idata%iostat
            info%stat = keep%idata%stat
            info%flag = flag
            if (flag == -17) info%flag = -16
            go to 25
          end if
          keep%used = keep%used - size(keep%imain)
          deallocate (keep%imain,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Try switching reals to a file
        if (keep%index(2) < 0) then
          keep%index(2) = -keep%index(2)
          loc1 = 1
          lw = int(keep%posfac,short) - 1
          call of01_write(keep%index(2),loc1,lw,keep%rmain, &
             flag,keep%rdata,lp=-1)
          if (flag < 0) then
            info%iostat = keep%rdata%iostat
            info%stat = keep%rdata%stat
            info%flag = flag
            if (flag == -17) info%flag = -16
            go to 25
          end if
          keep%used = keep%used - size(keep%rmain)*2
          deallocate (keep%rmain,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! In the indefinite case, nothing
! yet written to the in-core arrays so switch to files
        if (.not. keep%pos_def) then
          if (keep%index(4) < 0) then
            deallocate(keep%rwdelay,stat=st)
            keep%index(4) = -keep%index(4)
            keep%used = keep%used - 2*storage_indef
          end if
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Now try to write stack to file (at moment, nothing in the stack)
        if (keep%index(3) < 0) then
          keep%index(3) = -keep%index(3)
          keep%used = keep%used - size(keep%rwork)*2
          deallocate (keep%rwork,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Have not been able to allocate front
       info%flag = -30; info%stat = st
       call MA77_print_iflag(keep,nout,info%flag)
       go to 25
      end if
    end if
! frontal matrix has been allocated
 10 continue
    size_fa = lfa
!  write (6,*) 'keep%index ',keep%index(1:4)

! set keep%rfree and keep%ifree to point to where we will start writing
! factor data. Allow space for the entries of D to be held
    posdiag = keep%posfac
    inactive_fac = keep%posfac + n
    if (.not. keep%pos_def) inactive_fac = inactive_fac  + n
    keep%rfree = inactive_fac
    keep%ifree = keep%posint

! If arrays are being used in place of files, we have to ensure
! that the entries keep%posfac:keep%rfree contain meaningful data
! (otherwise, problem if we switch to using a file)
    if (keep%index(2) < 0) then
! use w54 as temporary array that we fill with zeros
      deallocate (w54,stat=st)
      lw54 = n
      if (.not. keep%pos_def) lw54 = lw54 + n
      allocate (w54(lw54),stat=st)
      if (st /= 0) go to 20
      w54(1:lw54) = zero
      call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain,       &
           posdiag,lw54,w54,flag,keep%rdata,-1,keep%maxstore,keep%used)
      if (flag < 0) go to 15
    end if

    deallocate (w54,stat=st)
    deallocate (w64,stat=st)
    if (keep%pos_def) then
      lw54 = max(nb,maxlen)*nb
      if (lsolve .and. nrhs >= 4) lw54 = lw54 + keep%maxelim*nb
      allocate (w54(lw54),stat=st)
      if (st /= 0) go to 20
    else
      lw64 = maxlen + maxlen*nb
      allocate (w64(lw64),stat=st)
      if (st /= 0) go to 20
    end if
! keep%maxelim_actual will hold actual max. number of eliminations at a node
! (only used in indef case but set in pos def case so not undefined)
    keep%maxelim_actual = 0

    deallocate (ip,stat=st)
    allocate (ip(maxlen),stat=st)
    if (st /= 0) go to 20

    if (.not.keep%pos_def) then
      deallocate (keep%size_ind,stat=st)
      deallocate (dpos,stat=st)
      deallocate (delay,stat=st)
      allocate (keep%size_ind(keep%l1:keep%l2),dpos(keep%l1:keep%l2), &
        delay(2,keep%l1:keep%l2),stat=st)
      if (st /= 0) go to 20
      keep%size_ind = 0
      dpos = 0_long
      delay = 0
    end if
    lreals = max(mvar,maxlen)

    deallocate (rlist,stat=st)
    deallocate (varlist,stat=st)
    deallocate (reals,stat=st)
    deallocate (pos,stat=st)
    deallocate (map,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)

    allocate (count(0:keep%maxdepth),cnode(1:keep%maxdepth),stat=st)
    if (st /= 0) goto 20

    allocate (reals(lreals),rlist(maxlen+1),varlist(maxlen),  &
              map(maxlen),pos(1:keep%n),stat=st)
    if (st /= 0) go to 20
! Initialise the array pos to zero
    pos = 0

    if (lsolve) then
      deallocate (xlocal,stat=st)
      allocate (xlocal(maxlen,nrhs),stat=st)
      if (st /= 0) go to 20
    end if

! initialise
    info%detlog  =  zero
    info%detsign  =  1
    info%ndelay = 0
    info%num_neg = 0
    info%num_nothresh = 0
    info%num_perturbed = 0
    info%matrix_rank = n - keep%null
    info%nfactor = 0_long
    info%nflops = 0_long
    info%ntwo = 0
    info%u = cntl64%u

!    info%wasted = 0

    keep%maxfa = 0_long
    keep%maxfrontb = 0
    keep%rtopmx = 1_long
    keep%rtopdmx = 1_long

    nroot = size(keep%roots)

    do ir = 1, nroot
      root = keep%roots(ir)
! Perform factorization for tree rooted at root.
! Initialise the top of the stack
      rtop = 1;  rtopd = 1
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. Assemble contributions
! from any children of node that are leaf nodes
! and then perform the partial factorization of node.
! Then either stack the contribution from node or
! assemble it into its parent (node is the ic-th child of its parent)
         if (pos_def) then
           if (present(scale)) then
             call factorize54(node,scale)
           else
             call factorize54(node)
           end if
           if (info%flag < 0 .or. depth <= 1) exit
           depth = depth - 1
           pnode = cnode(depth)
           ic = count(depth)
           call assemble54(pnode,node,ic,rtop)
         else
           if (present(scale)) then
             call factorize64(node,rtopd,scale)
           else
             call factorize64(node,rtopd)
           end if
           if (info%flag < 0 .or. depth <= 1) exit
           depth = depth - 1
           pnode = cnode(depth)
           ic = count(depth)
           call assemble64(pnode,node,ic,rtop)
         end if
         if (info%flag < 0) exit

! Go to parent
         node = pnode

      end do ! tree

! If pos. def. case indicated but found to be indefinite, allow code to
! rerun by not seeing keep%status negative.
      if (info%flag == -11 .and. pos_def) then
        call MA77_print_iflag(keep,nout,info%flag,st=info%stat, &
             ios=info%iostat)
        go to 25
      end if

      if (info%flag < 0) then
        keep%status = -1
        call MA77_print_iflag(keep,nout,info%flag,st=info%stat, &
             ios=info%iostat)
        go to 25
      end if
    end do ! roots

    keep%dfree = posdiag

    info%storage(1) = max(keep%mx_ifree,keep%ifree-1)
    info%storage(2) = keep%rfree - 1
    info%storage(3) = keep%rtopmx
    if (.not.pos_def) info%storage(4) = keep%rtopdmx
    info%minstore = info%storage(1) + &
      2*(info%storage(2) + info%storage(3) + info%storage(4))
    info%maxfront = keep%maxfrontb
    info%index(1:4) = keep%index(1:4)

    info%nio_read(1) = keep%idata%nio_read - info%nio_read(1)
    info%nio_read(2) = keep%rdata%nio_read - info%nio_read(2)
    info%nio_write(1) = keep%idata%nio_write - info%nio_write(1)
    info%nio_write(2) = keep%rdata%nio_write - info%nio_write(2)
    info%nwd_read(1) = keep%idata%nwd_read - info%nwd_read(1)
    info%nwd_read(2) = keep%rdata%nwd_read - info%nwd_read(2)
    info%nwd_write(1) = keep%idata%nwd_write - info%nwd_write(1)
    info%nwd_write(2) = keep%rdata%nwd_write - info%nwd_write(2)

    info%u = cntl64%u
    if (info%num_nothresh > 0) info%u = zero

    if (pos_def) then
! superfiles with index(4) have not been used
       info%index(4) = -keep%index(4)
    end if

    if (info%flag == 4) call MA77_print_iflag(keep,nout1,info%flag)

    if (info%flag >= 0) go to 25

15  info%iostat = keep%rdata%iostat
    info%stat = keep%rdata%stat
    info%flag = flag
    if (flag == -17) info%flag = -16
    go to 25

20  if (st /= 0) then
      info%flag = -1; keep%status = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
    end if

 25 deallocate (reals,stat=st)
    deallocate (map,stat=st)
    deallocate (pos,stat=st)
    deallocate (rlist,stat=st)
    deallocate (varlist,stat=st)
    deallocate (fa,stat=st)
    deallocate (xlocal,stat=st)
    deallocate (w64,stat=st)
    deallocate (dpos,stat=st)
    deallocate (delay,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (ip,stat=st)

    keep%flag = info%flag
    if (info%flag < 0) return
    keep%status = 3

! Do back subs if lsolve = .true.
    if (lsolve) then

      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
        write (control%unit_diagnostics,'(/a,i8,a)') &
      ' Factorisation complete. Solving for ',nrhs,' right-hand sides'
      temp(1:2) = info%nio_read(1:2)
      temp(3:4) = info%nio_write(1:2)
      temp(5:6) = info%nwd_read(1:2)
      temp(7:8) = info%nwd_write(1:2)

      if (present(scale)) then
        if (pos_def) then
          call MA77_solve(nrhs,lx,x,keep,control,info,scale=scale,job=3)
        else
          call MA77_solve(nrhs,lx,x,keep,control,info,scale=scale,job=4)
        end if
      else
        if (pos_def) then
          call MA77_solve(nrhs,lx,x,keep,control,info,job=3)
        else
          call MA77_solve(nrhs,lx,x,keep,control,info,job=4)
        end if
      end if
      deallocate (w54,stat=st)

      info%nio_read(1:2)  = temp(1:2) + info%nio_read(1:2)
      info%nio_write(1:2) = temp(3:4) + info%nio_write(1:2)
      info%nwd_read(1:2)  = temp(5:6) + info%nwd_read(1:2)
      info%nwd_write(1:2) = temp(7:8) + info%nwd_write(1:2)
      if (info%flag == 0) info%flag = keep%flag
    end if

!    write (6,'(A,F12.3)') 'CPU time for calls to ma54_factor  = ', tfactor
!    write (6,'(A,F12.3)') 'Wall time for calls to ma54_factor = ', stfactor

!    write (6,'(A,F12.3)') 'CPU time for calls to block form  = ', t54
!    write (6,'(A,F12.3)') 'Wall time for calls to block form = ', st54

    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      if (lsolve) then
        write (control%unit_diagnostics,'(/a)') &
       ' Completed factorisation and solve with:'
      else
        write (control%unit_diagnostics,'(/a)') &
       ' Completed factorisation with:'
      end if
      write (control%unit_diagnostics, &
     '(a,2(/a,i12),3(/a,es12.4),4(/a/2es12.4),5(/a,i12))') &
     ' information parameters (info%) :', &
     ' flag                Error flag                               = ', &
       info%flag, &
     ' maxfront            Maximum frontsize                        = ', &
       info%maxfront, &
     ' nfactor             Number of entries in L                   = ', &
       real(info%nfactor), &
     ' nflops              Number of flops performed                = ', &
       real(info%nflops), &
     ' minstore            Storage used in superfiles               = ', &
       real(info%minstore), &
 ' nio_read(1:2)       Number of records read from disk by OF01_read   = ',&
       real(info%nio_read(1:2)), &
 ' nio_write(1:2)      Number of records written to disk by OF01_write = ',&
       real(info%nio_write(1:2)), &
     ' nwd_read(1:2)       Number of scalars read by OF01_read      = ', &
       real(info%nwd_read(1:2)), &
     ' nwd_write(1:2)      Number of scalars written by OF01_write  = ', &
       real(info%nwd_write(1:2)), &
     ' ntwo                Number of 2x2 pivots used                = ', &
       info%ntwo, &
     ' ndelay              Number of delayed eliminations           = ', &
       info%ndelay, &
     ' rank                Computed rank                            = ', &
       info%matrix_rank, &
     ' num_neg             Computed number of negative eigenvalues  = ', &
       info%num_neg, &
     ' num_nothresh Pivots did not satisfy input threshold criteria = ', &
       info%num_nothresh

       if (info%num_perturbed == 0) then
         write (control%unit_diagnostics,'(a,es12.4)') &
     ' threshold parameter used (info%u)                       = ', &
         info%u
       else
          write (control%unit_diagnostics,'(a,i12)') &
     ' number of perturbed pivots (info%num_perturbed)              = ', &
         info%num_perturbed
       end if
    end if

contains

!*******************************

   subroutine factorize54(node,scale)
! This subroutine assembles contributions from any children of node
! that are leaf nodes and then performs the partial factorization at node
! (pos definite case).

    integer(short), intent(in) :: node ! node in tree
    real(wp), intent(in), optional :: scale(*)

! Possible error returns:
!  -1   Allocation error
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -11   Problem found to be not pos.def.
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -22   For one or more elements/rows, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.

      nvar = abs(keep%size(node))
      if (nvar == 0) return
! If node is a leaf, immediate return.
      if (node <= keep%nelt) return

      nvar_long = nvar
      nfront = nvar;  nfront_long = nvar
      ln = (nfront_long*(nfront_long+1))/2
! set ip1 to point to first location in fa that will be used for front at node
      ip1 = 1 + size_fa - ln
      lfa = ln

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild
! In element case, add up the number of children of node that are leaf nodes
! .... they are the last children (in row case, leaf nodes are original
! rows and they are not part of the tree)
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if
! If all the children of node are leaf children, initialize front
      if (non_leaf == 0) fa(ip1:ip1+lfa-1) = zero

! Read the variable list for node into rlist.
      loc = keep%ifile(node)
      call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
           rlist,flag,keep%idata,-1)
      if (flag < 0) go to 480
      do i = 1,nvar
        ivar = rlist(i)
        pos(ivar) = i
      end do

! In element case, must now deal with children that are leaf nodes
! and in row case, deal with original rows

! Set ip(j) to point to start of col j in array fa.
      ip(1) = ip1
      do j = 2,nvar
        ip(j) = ip(j-1) + nvar - j + 2
      end do

      if (keep%element_input) then
        do ic = non_leaf+1,nchild
          jcnode = keep%tree(node)%child(ic)
! Check reals were entered.
          if (keep%rfile(jcnode) == -1 .and. keep%size(jcnode) /= 0) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
          cnvar = abs(keep%size(jcnode))
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
               varlist,flag,keep%idata,-1)
          if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos.
          do j = 1,cnvar
            k = varlist(j)
            map(j) = pos(k)
          end do

! Read reals from main file one column at a time
          loc = cnvar
          loc = keep%rfile(jcnode) - (loc*loc+loc)/2
          do j = 1,cnvar
            lreals = cnvar - j + 1
            call MA77_read_real(keep%index(2),keep%rmain,loc,lreals, &
                 reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480

! Add the reals from the original element into the front (optionally scale)
            k = 1
            jj = map(j)
            if (present(scale)) then
              k1 = varlist(j)
              do i = j, cnvar
                ii = map(i)
                if (ii >= jj) then
                  klong = ip(jj) + ii - jj
                else
                  klong = ip(ii) + jj - ii
                end if
                k2 = varlist(i)
                reals(k) = scale(k1)*reals(k)*scale(k2)
                fa(klong) = fa(klong) + reals(k)
                k = k + 1
              end do
            else
              do i = j, cnvar
                ii = map(i)
                if (ii >= jj) then
                  klong = ip(jj) + ii - jj
                else
                  klong = ip(ii) + jj - ii
                end if
                fa(klong) = fa(klong) + reals(k)
                k = k + 1
              end do
            end if
            loc = loc + lreals
          end do
        end do

      else
! For row entry, add into the front the original rows corresponding
! to the eliminated variables.
        do j = 1, keep%tree(node)%nelim
          jcnode = rlist(j)
          cnvar = abs(keep%size(jcnode))
          if (cnvar == 0) cycle
          if (keep%rfile(jcnode) == -1) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
! Read in the variable list for row jcnode and its reals
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar,varlist,&
               flag,keep%idata,-1)
          if (flag < 0) go to 480
! Recall that keep%rfile(jcnode) points to position in the main reals
! file after the last entry in the list of reals for jcnode
          loc = keep%rfile(jcnode) - cnvar
          call MA77_read_real(keep%index(2),keep%rmain,loc,cnvar,reals,flag,&
               keep%rdata,-1)
          if (flag < 0) go to 480

! Loop over entries in the row jcnode, putting
! them into row/col j of the frontal matrix (lower triangular part).
! ipp points to start of col. j
           jlong = j
           ipp = ip1 + (nvar_long+1)*(jlong-1) - (jlong*(jlong-1))/2
           if (present(scale)) then
             k1 = jcnode
             do jj = 1, cnvar
               k2 = varlist(jj)
               i = pos(k2)
               if (i < j) cycle
! Put entry into col j, row i.
               reals(jj) = scale(k1)*reals(jj)*scale(k2)
               fa(ipp+i-j) = fa(ipp+i-j) + reals(jj)
             end do
           else
             do jj = 1, cnvar
               k2 = varlist(jj)
               i = pos(k2)
               if (i < j) cycle
! Put entry into col j, row i.
               fa(ipp+i-j) = fa(ipp+i-j) + reals(jj)
             end do
           end if
        end do
      end if

! Restore pos to zero
      do i = 1,nvar
        ivar = rlist(i)
        pos(ivar) = 0
      end do

! Have now assembled the frontal matrix and are ready to factorize it.
! It is of order nfront=nvar, held in lower triangular packed form

     nelim = keep%tree(node)%nelim
     keep%maxfrontb = max(keep%maxfrontb,nfront)

!     call cpu_time(t1)
!     call mysystem_clock(st1)

! Rearrange to block hybrid form and Cholesky factorization of nelim columns

     call ma54_to_block(nfront,nelim,nb,fa(ip1:ip1+lfa-1),w54,flag)

     call ma54_factor(nfront,nelim,nb,fa(ip1:ip1+lfa-1),w54,flag)

     if (flag > 0) then
! Matrix not positive definite. error return -11
!      write (6,*) 'error from ma54_factor. flag = ',flag,nfront,nelim
       info%flag = -11
       go to 500
     end if

     lreals = (nelim*nelim+nelim)/2 + nelim*(nfront-nelim)
!      write (6,'(a,6i6)')  &
!     'node,nfront,ncand,nelim,delay ',node,nfront,nelim,nelim,nelim-nelim

     if (nelim > 0) then

! obtain the diagonal entries and store in main real file
       call ma54_diag(nfront,nelim,nb,fa(ip1:ip1+lfa-1),w54)
       locw = posdiag
       call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain, &
            locw,nelim,w54,flag,keep%rdata,-1,keep%maxstore, &
            keep%used,inactive=keep%posfac)
       if (flag < 0) go to 480
       posdiag = posdiag + nelim
! Compute the determinant
       do k = 1,nelim
         info%detlog = info%detlog + log(abs(w54(k)))
         if (w54(k) < 0) info%detsign = -info%detsign
       end do

! Store factor in fa(ip1:ip1+lreals-1) to the main file.
       lreals = nfront
       ipp = ip1
       locw = keep%rfree
       do j = 1,nelim
         call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain, &
              locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1, &
              keep%maxstore,keep%used,inactive=inactive_fac)
         if (flag < 0) go to 480
         locw = locw + lreals
         ipp = ipp + lreals
         lreals = lreals - 1
       end do

! Update keep%maxfa to hold the longest factor entry written
       llreals = nelim
       llreals = (llreals*llreals+llreals)/2 + llreals*(nfront-llreals)
       keep%maxfa = max(keep%maxfa,llreals)

! Set info%nfactor to hold no. of entries in factor (long integers needed)
       lnj = nelim
       lnj1 = nfront - nelim
       llreals = (lnj*lnj+lnj)/2_long + lnj*lnj1
       keep%rfree = keep%rfree + llreals
       info%nfactor = info%nfactor + llreals

      lnj = nfront
      info%nflops = info%nflops + &
            (lnj*(lnj+1)*(2_long*lnj+1))/6_long
      lnj = nfront - nelim
      info%nflops = info%nflops - &
            (lnj*(lnj+1)*(2_long*lnj+1))/6_long

     end if

! Partial factorization complete.

!!!!!!!
! If x was supplied then perform the forward subs.
     if (lsolve .and. nelim > 0) then
       if (size(xlocal,1) < nfront) then
         deallocate (xlocal,stat=st)
         allocate (xlocal(nfront,nrhs),stat=st)
         if (st /= 0) go to 490
       end if
       do i = 1,nfront
         jvar = rlist(i)
         xlocal(i,1:nrhs) = x(jvar,1:nrhs)
       end do
       if (nrhs >= 4) then
         call ma54_forward2(nfront,nelim,nb,nrhs,fa(ip1:ip1+lfa-1), &
              xlocal,maxlen,nb,w54,flag)
       else
         do j = 1, nrhs
           call ma54_forward1(nfront,nelim,nb,fa(ip1:ip1+lfa-1),    &
                xlocal(:,j),flag)
         end do
       end if
       do i = 1, nfront
         jvar = rlist(i)
         x(jvar,1:nrhs) = xlocal(i,1:nrhs)
       end do
!       if (flag /= 0) then
!         write (6,*) 'error forward: flag=', flag
!         stop
!       end if
     end if

!!!!!
      if (info%flag >= 0) go to 500

480   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16 
      go to 500

490   info%flag = -1
      info%stat = st
      go to 500

500   continue
      info%index(1:4) = keep%index(1:4)
! files not used for work array
      info%index(4) = -keep%index(4)
      return

    end subroutine factorize54
!*******************************

    subroutine assemble54(node,cnode,ic,rtop)

! This subroutine either stacks the contribution from cnode
! or assembles it into its parent node. cnode is the ic-th child of node.

    integer(short), intent(in) :: node ! node in tree
    integer(short), intent(in) :: cnode ! child of node
    integer(short), intent(in) :: ic ! position of ic in sibling list
    integer(long), intent(inout) :: rtop ! points to top of main real stack

! Possible error returns:
!  -1   Allocation error
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild
! In element case, add up the number of children of node that are leaf nodes
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if

! set ipc to point to first location in fa that was used for cnode
        cnvar = keep%size(cnode) - keep%tree(cnode)%nelim
        length = cnvar
        length = (length*length + length)/2
        ipc = 1 + size_fa - length

        nvar = abs(keep%size(node))
        nfront = nvar
        nvar_long = nvar
        ln = (nvar_long*nvar_long + nvar_long)/2
! set ip1 to point to first location in fa that will be used for front at node
        ip1 = 1 + size_fa - ln
        lfa = ln
        splitp = keep%splitp(node)
! If we are not yet at the split point, write generated element for cnode onto
! top of the real stack and return
        if (ic < splitp) then
          locw = rtop
          lreals = cnvar
          ipp = ipc
          do j = 1,cnvar
            call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork, &
                 locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,        &
                 keep%maxstore,keep%used)
            if (flag < 0) go to 480
            ipp = ipp + lreals
            locw = locw + lreals
            lreals = lreals - 1
          end do
! Adjust top of stack
          rtop = rtop + length
          keep%rtopmx = max(keep%rtopmx,rtop)
          return
        end if

! Read the variables for node from the main integer superfile.
        loc = keep%ifile(node)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
             rlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
        do i = 1,nvar
          ivar = rlist(i)
          pos(ivar) = i
        end do

! Set ip(j) to point to start of col j in array fa.
        ip(1) = ip1
        do j = 2,nvar
          ip(j) = ip(j-1) + nvar - j + 2
        end do

! read integer list for cnode from main integer file
        loc = keep%ifile(cnode) + keep%tree(cnode)%nelim
        call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
             varlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
! Set mapping from cnode into the front
        do j = 1,cnvar
          k = varlist(j)
          map(j) = pos(k)
        end do

        if (ic == splitp .or. ic == non_leaf) then
          if (cnvar /= nvar) then
! Expand cnode into the front
            j1 = 1
            do j = 1,cnvar
              j2 = map(j)
! col. j of cnode maps to col. j2 of front
! set cols j1:j2-1 of front to zero
              fa(ip(j1):ip(j2)-1) = zero
! Expand col. j of cnode into col. j2 of front.
! col j of cnode is currently at col k=j+nvar-cnvar of front
              k = j + nvar - cnvar
              if (j2 /= k) then
                ipj2 = ip(j2)
                ipk = ip(k)
                fa(ipj2:ip(j2+1)-1) = zero
                do l = j,cnvar
                  i = map(l)
                  fa(ipj2+i-j2) = fa(ipk+l-j)
                end do
              else
! j2 = k, so all columns are now in place
                j1 = nvar + 1
                exit
              end if
              j1 = j2 + 1
            end do
! Check final cols are set to 0.
            klong = size_fa
            if (j1 <= nvar) fa(ip(j1):klong) = zero
          end if
        end if

        if (ic == splitp) then
! At split point.
! Merge children 1:splitp-1 into rows/columns 1:nfront of frontal matrix
          do 160 jc = splitp-1, 1, -1
            jcnode = keep%tree(node)%child(jc)
            nelim_jc = keep%tree(jcnode)%nelim
            loc = keep%ifile(jcnode) + nelim_jc
            cnvar = keep%size(jcnode) - nelim_jc
            call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
                 varlist,flag,keep%idata,-1)
            if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos (length n)
            do j = 1,cnvar
              k = varlist(j)
              map(j) = pos(k)
            end do
! Read in reals for child from the top of stack, one col. at a time,
! using map to read directly into the front.
            length = cnvar
            length = (length*length + length)/2
            loc = rtop - length
            lreals = cnvar
! move pointer to top of stack
            rtop = loc
            do j = 1, cnvar
              ilong = map(j)
              ipp = ip1 + nvar_long*(ilong-1) - (ilong*(ilong-1))/2
              call MA77_read_discard_real(keep%index(3),keep%rwork,loc,lreals,&
                   fa(ipp:ipp+nvar-1),flag,keep%rdata,-1,map=map(j:cnvar))
              if (flag < 0) go to 480
              loc = loc + lreals
              lreals = lreals - 1
            end do
! End of loop over children
160       continue

! At this point fa holds the reals in the frontal matrix in packed form.
! If split point is not also last non-leaf child, write fa onto the stack
          if (ic /= non_leaf) then
            lreals = nvar
            loc = rtop
            ipp = ip1
            do j = 1,nvar
              call MA77_write_real(keep%index(3),keep%size_rmain,keep%rwork, &
                   loc,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,         &
                   keep%maxstore,keep%used)
              if (flag < 0) go to 480
              ipp = ipp + lreals
              loc = loc + lreals
              lreals = lreals - 1
            end do
            rtop = rtop + (nvar_long*nvar_long+nvar_long)/2
            keep%rtopmx = max(keep%rtopmx,rtop)
          end if

        else if (ic < non_leaf) then
! We have passed split point and not yet reached last non-leaf child.
! The reals for jcnode are still in fa.
! Loop over the cols of cnode, merging the child into its parent.
            do j = 1,cnvar
              i = map(j)
! Read column i of the frontal matrix from the stack. Find its position.
               loc = rtop - (nvar*nvar+nvar)/2
! loc now points to the first entry in the front on the stack.
! Set loc to point to first entry in col. i and
! lreals to number of entries in col. i.
              loc = loc + nvar*(i-1) - ((i-2)*(i-1))/2
              lreals = nvar - i + 1
              call MA77_read_real(keep%index(3),keep%rwork,loc,lreals,reals,&
                   flag,keep%rdata,-1)
              if (flag < 0) go to 480
! ipp points to start of col. j of jcnode, which is held in fa
              jlong = j
              ipp = ipc + cnvar*(jlong-1) - ((jlong-2)*(jlong-1))/2
              do k = j,cnvar
                kk = map(k)
                reals(1+kk-i) = reals(1+kk-i) + fa(ipp+k-j)
              end do
! Write the updated reals back onto the stack
              call MA77_write_real(keep%index(3),keep%size_rmain,keep%rwork, &
                  loc,lreals,reals,flag,keep%rdata,-1,keep%maxstore,keep%used)
              if (flag < 0) go to 480
            end do

        else if (ic == non_leaf) then
! Dealing with last non-leaf child (and it is not split point)
! Read frontal matrix from top of stack, one column at a time and add
! in with last non-leaf child (which we have already expanded in fa).
          loc = rtop - (nvar_long*nvar_long+nvar_long)/2
          rtop = loc
          ipp = ip1
          lreals = nvar
          do j = 1,nvar
            call MA77_read_discard_real(keep%index(3),keep%rwork,loc,&
                 lreals,reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480
            fa(ipp:ipp+lreals-1) = fa(ipp:ipp+lreals-1) + reals(1:lreals)
            ipp = ipp + lreals
            loc = loc + lreals
            lreals = lreals - 1
          end do

        end if

      if (info%flag >= 0) return

480   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16 

    end subroutine assemble54

!***********************************************

   subroutine factorize64(node,rtopd,scale)
! This subroutine assembles contributions from any children of node
! that are leaf nodes, adds in contributions from delayed
! pivots and then performs the partial factorization at node
! (indefinite case).

    integer(short), intent(in) :: node ! node in tree
    integer(long) :: rtopd ! points to the top of the real stack for
!                            delayed pivots (only needed in indefinite case)
    real(wp), intent(in),optional :: scale(*)

! Possible error returns:
!  -1   Allocation error
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -11   Problem found to be singular.
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -22   For one or more elements/rows, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.
! -30   Front size too large

      nvar = abs(keep%size(node))
      if (nvar == 0) return
! If node is a leaf, immediate return.
      if (node <= keep%nelt) return

      nvar_long = nvar
      ln = (nvar_long*(nvar_long+1))/2
      ip1 = 1 + size_fa - ln

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild

! In element case, add up the number of children that are leaf nodes
! .... they are the last children (in row case, leaf nodes are original
! rows and they are not part of the tree)
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if

! If node has no non-leaf children, initialize front
     cdelay = 0
     if (non_leaf == 0) then
       fa(ip1:ip1+ln-1) = zero
     else
! Add up number of delayed cols being passed up the children
       do j = 1,non_leaf
         jcnode = keep%tree(node)%child(j)
         cdelay = cdelay + delay(2,jcnode)
       end do
     end if
     nfront = nvar + cdelay
     nfront_long = nfront
     lfa = (nfront_long*(nfront_long+nb+1))/2

! Read the variables for node from the main integer superfile into
! the end of rlist and set mapping into front
      loc = keep%ifile(node)
      call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
           rlist(cdelay+1:cdelay+nvar),flag,keep%idata,-1)
      if (flag < 0) go to 480
      do i = 1,nvar
        ivar = abs(rlist(cdelay+i))
        pos(ivar) = i
      end do

! In element case, must now deal with children that are leaf nodes
! and in row case, deal with original rows

! Set ip(j) to point to start of jth col.
      ip(1) = ip1
      do j = 2,nvar
        ip(j) = ip(j-1) + nvar - j + 2
      end do
      if (keep%element_input) then
        do ic = non_leaf+1,nchild
          jcnode = keep%tree(node)%child(ic)
! Check reals were entered.
          if (keep%rfile(jcnode) == -1 .and. keep%size(jcnode) /= 0) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
          cnvar = abs(keep%size(jcnode))
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
               varlist,flag,keep%idata,-1)
          if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos.
          do j = 1,cnvar
            k = varlist(j)
            map(j) = pos(k)
          end do

! Read reals from main file one column at a time
          loc = cnvar
          loc = keep%rfile(jcnode) - (loc*loc+loc)/2
          do j = 1,cnvar
            inelrs = cnvar - j + 1
            call MA77_read_real(keep%index(2),keep%rmain,loc,inelrs, &
                 reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480

! Add the reals from the original element into the front
            k = 1
            jj = map(j)
            if (present(scale)) then
              k1 = varlist(j)
              do i = j, cnvar
                ii = map(i)
                if (ii >= jj) then
                  klong = ip(jj) + ii - jj
                else
                  klong = ip(ii) + jj - ii
                end if
                k2 = varlist(i)
                reals(k) = scale(k1)*reals(k)*scale(k2)
                fa(klong) = fa(klong) + reals(k)
                k = k + 1
              end do
            else
              do i = j, cnvar
                ii = map(i)
                if (ii >= jj) then
                  klong = ip(jj) + ii - jj
                else
                  klong = ip(ii) + jj - ii
                end if
                fa(klong) = fa(klong) + reals(k)
                k = k + 1
              end do
            end if
            loc = loc + inelrs
          end do
        end do

      else
! For row entry, add into the front the original rows corresponding
! to the eliminated variables.
        do j = 1, keep%tree(node)%nelim
          jcnode = abs(rlist(cdelay+j))
          cnvar = abs(keep%size(jcnode))
          if (cnvar == 0) cycle
          if (keep%rfile(jcnode) == -1) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
! Read in the variable list for row jcnode and its reals
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar,varlist,&
               flag,keep%idata,-1)
          if (flag < 0) go to 480
! Recall that keep%rfile(jcnode) points to position in the main reals
! file after the last entry in the list of reals for jcnode
          loc = keep%rfile(jcnode) - cnvar
          call MA77_read_real(keep%index(2),keep%rmain,loc,cnvar,reals,flag,&
               keep%rdata,-1)
          if (flag < 0) go to 480

! Loop over entries in the row jcnode, putting
! them into col j of the frontal matrix (lower triangular part).
! ipp points to start of col. j
           jlong = j
           ipp = ip1 + (nvar_long+1)*(jlong-1) - (jlong*(jlong-1))/2
           if (present(scale)) then
             k1 = jcnode
             do jj = 1, cnvar
               k2 = varlist(jj)
               i = pos(k2)
               if (i < j) cycle
! Put entry into col j, row i.
               reals(jj) = scale(k1)*reals(jj)*scale(k2)
               fa(ipp+i-j) = fa(ipp+i-j) + reals(jj)
             end do
           else
             do jj = 1, cnvar
               k2 = varlist(jj)
               i = pos(k2)
               if (i < j) cycle
! Put entry into col j, row i.
               fa(ipp+i-j) = fa(ipp+i-j) + reals(jj)
             end do
           end if
        end do
      end if

! All children dealt with.
! We have to put the delayed columns into first columns of frontal matrix.
! Read in integer and real data corresponding to delayed pivots
! from the integer and real delay stacks
      nchild = size(keep%tree(node)%child)
      ld = 1
      nfront_long = nfront
      ln = (nfront_long*(nfront_long+1))/2
      ip1 = 1 + size_fa - ln
      ipp = ip1
      k = 1
!     write (6,*) 'node,nfront,size_fa,ln,ip1', &
!                  node,nfront,size_fa,ln,ip1
      do jc = nchild, 1, -1
        jcnode = keep%tree(node)%child(jc)
        if (jcnode <= keep%nelt) cycle
        ldc = delay(2,jcnode)
! cycle if no delayed pivots for jcnode.
        if (ldc == 0) cycle
        cnvar = delay(1,jcnode)

! Read delayed variables from integer file into varlist
        loci = dpos(jcnode)
        call MA77_read_integer(keep%index(1),keep%imain,loci,cnvar,varlist,&
             flag,keep%idata,-1)
        if (flag < 0) go to 480
!  write (6,*) 'read: ldc,loci,varlist',ldc,loci,varlist(1:cnvar)
! Copy indices for the delayed pivots at jcnode into start of rlist
        rlist(k:k+ldc-1) = varlist(1:ldc)
        k = k + ldc
! Set mapping between the other variables in delayed cols and current front

        if (cnvar > size(map)) then
          deallocate (map,stat=st)
          allocate (map(cnvar),stat=st)
          if (st /= 0) go to 490
        end if
        do j = ldc+1,cnvar
          jj = varlist(j)
          map(j) = pos(jj)
        end do
        lreals = cnvar
        ilong = ldc
        locr = rtopd - (ilong*(cnvar-ilong) + (ilong*ilong+ilong)/2)
! Reset top of real stack and then read in reals one column at a time
! and put into the start of the frontal matrix.
        rtopd = locr
        do i = 1, ldc
          call MA77_read_discard_real(keep%index(4),keep%rwdelay,locr,&
               lreals,reals,flag,keep%rdata,-1)
          if (flag < 0) go to 480
! ipp is start of column ld in array fa
!         write (6,*) 'ipp,ldc,i,ld,ipp+ldc-i+1,ipp+nfront-ld', &
!                      ipp,ldc,i,ld,ipp+ldc-i+1,ipp+nfront-ld
          fa(ipp:ipp+ldc-i) = reals(1:ldc-i+1)
! Set zeros in the rest of the col. to zero
          fa(ipp+ldc-i+1:ipp+nfront-ld) = zero

! Use map to put the remaining entries into the front
          ipp = ipp + cdelay - ld
          l = ldc - i + 1
          do ii = ldc+1,cnvar
            jrow = map(ii)
            l = l + 1
            fa(ipp+jrow) = reals(l)
          end do
          locr = locr + lreals
          lreals = lreals - 1
          ipp = ipp + nvar + 1
          ld = ld + 1
        end do

      end do

! Reset pos to zero
      do i = 1,nvar
        ivar = abs(rlist(cdelay+i))
        pos(ivar) = 0
      end do
!!!!(we want to use negative sign in rlist to indicate 2x2 pivot)

      continue

      keep%maxfrontb = max(keep%maxfrontb,nfront)

! Have now added the leaf nodes and the delayed
! columns into the frontal matrix and are ready to factorize it.
! It is of order nfront, held in lower triangular packed form
! Let ncand be no. of potential pivots

      ncand = keep%tree(node)%nelim + cdelay

! nelim will be the number of eliminations actually performed.
! reals is used to hold the inverse of the block diagonal
! in the LDL^T factorization. It must be length 2*ncand
      i = 2*ncand
      if (size(reals) < i) then
        deallocate (reals,stat=st)
        allocate(reals(max(i,int(multiplier*i))),stat=st)
        if (st /= 0) go to 490
      end if
      i = nfront + nfront*nb
      if (size(w64) < i) then
        deallocate (w64,stat=st)
        allocate(w64(max(i,int(multiplier*i))),stat=st)
        if (st /= 0) go to 490
      end if

! varlist will be used to hold the ma64 array perm
! If we want to flag up 2x2 pivots, then we do this using varlist here
! check whether 2x2 pivots were supplied
       if (keep%ntwo > 0) then
         do i = 1,cdelay
           varlist(i) = 1
         end do
         k = cdelay + 1
         do i = cdelay+1,ncand
           varlist(k) = rlist(k)
           if (rlist(k) < 0) then
! reset sign
             rlist(k) = -rlist(k)
! ensure negative flag for the partner
             varlist(k+1) = -rlist(k+1)
             k = k + 1
           end if
           k = k + 1
           if (k > ncand) exit
         end do
       end if

!      write (6,*) 'Before MA64: nfront,nvar,ncand,nb,nbi,lfa,cdelay',&
!                   nfront,nvar,ncand,nb,nbi,lfa,cdelay
!      write (6,*) 'rlist',rlist(1:nfront)
!      write (6,*) 'varlist',varlist(1:ncand)

!      write (6,*) 'ip1,ln,ip1+ln-lfa,size_fa',&
!                   ip1,ln,ip1+ln-lfa,size_fa
!      if (cdelay > 0) then
!       write (6,'(4(es12.4,2x))') fa(ip1:ip1+ln-1)
!       stop
!      end if

      call ma64_factor(nfront,ncand,nb,nbi,fa(ip1+ln-lfa:ip1+ln-1),lfa,  &
           cntl64,nelim,lq,varlist,reals,w64,info64,s=cdelay)

      if (info64%flag == -13) then
        info%flag = -40;  goto 500
      else if (info64%flag < 0) then
!        write (6,*) 'unexpected error from ma64_factor. flag = ',info64%flag
        info%flag = -99;  go to 500
      end if
      if (ncand == nfront) then
        if (nelim /= ncand) then
!         write (6,*) 'error from ma64_factor. not enough pivots chosen'
!         write (6,*) 'nelim,ncand,num_zero', nelim,ncand,info64%num_zero
         info%flag = -99; go to 500
        end if
      end if
!      if (ncand /= nelim) write (6,*) 'node,ncand,nelim',node,ncand,nelim
!      write (6,*) 'After  MA64: ncand,nelim',ncand,nelim

! Keep a record of the actual max. number of eliminations at a node
      keep%maxelim_actual = max(keep%maxelim_actual,nelim)

      if (info64%num_zero > 0) then
        info%matrix_rank = info%matrix_rank - info64%num_zero
      !  write (6,*) rank, num_zero',info%matrix_rank, info64%num_zero
        if (.not. control%action) then
          info%flag = -11;   go to 500
        else
          info%flag = 4
        end if
      end if

      info%num_nothresh  = info%num_nothresh + info64%num_nothresh
      info%num_perturbed = info%num_perturbed + info64%num_perturbed
      info%ntwo          = info%ntwo + info64%num_2x2
      info%num_neg       = info%num_neg + info64%num_neg
      info%detlog        = info%detlog + info64%detlog
      info%detsign       = info%detsign*info64%detsign

! If we have used a smaller u, we will use it at the next stage
      cntl64%u = info64%u

! Store factor held in fa(ip1+ln-lfa:ip1+ln-lfa+lq-1) to the main file
      if (nelim > 0) then
        ipp = 1 + size_fa - lfa   ! = ip1 + ln - lfa
        locw = keep%rfree
        lreals = nfront
        do j = 1,nelim
          call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain, &
               locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,        &
               keep%maxstore,keep%used,inactive=inactive_fac)
          if (flag < 0) go to 480
          locw = locw + lreals
          ipp = ipp + lreals
          info%nfactor = info%nfactor + lreals
          lreals = lreals - 1
        end do
        keep%rfree = locw

! set llreals to be the number of entries stored
        llreals = nelim
        llreals = (llreals*llreals+llreals)/2 + llreals*(nfront-llreals)

! Update keep%maxfa to hold upper bound on
! the longest factor entry written (we do not
! include the diagonal entries since we will read these later with
! separate array)
        keep%maxfa = max(keep%maxfa,llreals)

! Also store the block diagonal (held in reals)
        dreals = 2*nelim
        locw = posdiag
        call MA77_write_real(keep%index(2),keep%size_rmain,keep%rmain,       &
               locw,dreals,reals,flag,keep%rdata,-1,keep%maxstore,keep%used, &
               inactive=keep%posfac)
        if (flag < 0) go to 480
! 10 sept. 2008. No longer include D in the count info%nfactor
!        info%nfactor = info%nfactor + dreals
        posdiag = posdiag + dreals

! Flop count
        do j = 0,nelim-1
          lnj = nfront - j
          info%nflops = info%nflops + lnj*lnj
        end do
      end if

! Set delay(1,node) to hold frontsize after eliminations and
! delay(2,node) to hold number of delayed pivots
      delay(1,node) = nfront - nelim
      ldc = ncand - nelim
      delay(2,node) = ldc

! Write delayed pivots to the delayed integer
! and real stacks, and write integer data to main integer file

! Apply permutation computed by factorization to fully summed cols.
!       write (6,*) 'nelim,delay(1:2)',nelim,delay(1:2)
!       write (6,*) 'perm ',varlist(1:ncand)

! Use map to take temporary copy (because f95 compiler with -O is giving me
! wrong answer if I use following)
!       rlist(1:ncand) = rlist(varlist(1:ncand))
        if (size(map) < nfront) then
          deallocate (map,stat=st)
          allocate (map(nfront),stat=st)
          if (st /= 0) go to 490
        end if
        map(1:nfront) = rlist(1:nfront)
        do i = 1,ncand
          rlist(i) = map(varlist(i))
        end do

! Store no. of eliminations performed
        rlist(nfront+1) = nelim
! Write rlist onto the end of the main integer file
        locw = keep%ifree
        dpos(node) = locw + nelim
!       write (6,*) 'rlist,locw',rlist(1:nfront+1),locw
        call MA77_write_integer(keep%index(1),keep%size_imain,keep%imain, &
             locw,nfront+1,rlist,flag,keep%idata,-1,keep%maxstore,keep%used)
        if (flag < 0) go to 480
        locw = locw + nfront + 1
        keep%ifree = locw
        keep%size_ind(node) = nfront
        info%ndelay = info%ndelay + ldc

        if (ldc > 0) then
! Write the reals for delayed cols onto the temporary real stack
          nfront_long = nfront
          ipp = 1 + size_fa - ((nfront_long-nelim)*(nfront_long-nelim+1))/2
          locw = rtopd
          lreals = nfront - nelim
          ilong = ldc
          length = ilong*(nfront_long-ncand) + (ilong*ilong+ilong)/2
          do j = 1,ldc
            call MA77_write_real(keep%index(4),keep%size_rwdelay,keep%rwdelay,&
                 locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,  &
                 -1,keep%maxstore,keep%used)
            if (flag < 0) go to 480
            locw = locw + lreals
            ipp = ipp + lreals
            lreals = lreals - 1
          end do
          rtopd = rtopd + length
          keep%rtopdmx = max(keep%rtopdmx,rtopd)
        end if
!!!!!!!
! If x was supplied then perform the forward subs.
     if (lsolve .and. nelim > 0) then
       if (size(xlocal,1) < nfront) then
         deallocate (xlocal,stat=st)
         allocate (xlocal(nfront,nrhs),stat=st)
         if (st /= 0) go to 490
       end if
       do i = 1,nfront
         jvar = rlist(i)
         xlocal(i,1:nrhs) = x(jvar,1:nrhs)
       end do
       ipp = 1 + size_fa - lfa
       llreals = nelim
       llreals = (llreals*(2*nfront_long-llreals+nb))/2
       l = size(xlocal,1)
       if (nrhs >= 4) then
         call ma64_solveL2(nfront,nelim,nb,nrhs,xlocal,l, &
              flag,fa(ipp:ipp+llreals-1),lq)
       else
         do j = 1, nrhs
           call ma64_solveL1(nfront,nelim,nb,xlocal(:,j),flag, &
                fa(ipp:ipp+llreals-1),lq)
         end do
       end if

       do i = 1, nfront
         jvar = rlist(i)
         x(jvar,1:nrhs) = xlocal(i,1:nrhs)
       end do
     end if

!!!!!
      if (info%flag >= 0) go to 500

480   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16 
      go to 500

490   info%flag = -1
      info%stat = st
      go to 500

500   continue
      info%index(1:4) = keep%index(1:4)
      return

    end subroutine factorize64

!*******************************

    subroutine assemble64(node,cnode,ic,rtop)

! This subroutine either stacks the contribution from cnode
! or assembles it into its parent node. cnode is the ic-th child of node.

    integer(short), intent(in) :: node ! node in tree
    integer(short), intent(in) :: cnode ! child of node
    integer(short), intent(in) :: ic ! position of ic in sibling list
    integer(long), intent(inout) :: rtop ! points to top of main real stack

! Possible error returns:
!  -1   Allocation error
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -30   Front size too large

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild
! In element case, add up the number of children of node that are leaf nodes
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if

! set ipc to point to first location in fa that was used for cnode
        cnvar = keep%size(cnode) - keep%tree(cnode)%nelim
        length = cnvar
        length = (length*length+length)/2
        ipc = 1 + size_fa - length

        nvar = abs(keep%size(node))
        nfront = nvar
        nvar_long = nvar
! set ip1 to point to first location in fa that will be used for front at node
        lfa = (nvar_long*nvar_long + nvar_long)/2
        ip1 = 1 + size_fa - lfa 

! If we are not yet at the split point, write generated element for cnode onto
! top of the real stack and return
        splitp = keep%splitp(node)
        if (ic < splitp) then
          locw = rtop
          lreals = cnvar
          ipp = ipc
          do j = 1,cnvar
            call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork, &
                 locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,  &
                 keep%maxstore,keep%used)
            if (flag < 0) go to 480
            ipp = ipp + lreals
            locw = locw + lreals
            lreals = lreals - 1
          end do
! Adjust top of stack
          rtop = rtop + length
          keep%rtopmx = max(keep%rtopmx,rtop)
          return
        end if

        if (ic == non_leaf) then
! At the last non-leaf child, extra space may be needed in indefinite case
! to accommodate the delayed columns (they go at the start of
! the frontal matrix)
! Add up number of delayed cols and store max. number
          cdelay = 0
          d1max = 0
          do j = 1,non_leaf
            jcnode = keep%tree(node)%child(j)
            cdelay = cdelay + delay(2,jcnode)
            d1max = max(d1max,delay(1,jcnode))
          end do
          nfront = nvar + cdelay
          ncand = keep%tree(node)%nelim + cdelay
          if (size(rlist) < nfront+1) then
            deallocate (rlist,stat=st)
            allocate (rlist(int(nfront*multiplier)+1),stat=st)
            if (st /= 0) go to 490
          end if
          if (size(varlist) < max(d1max,ncand)) then
            deallocate (varlist,stat=st)
            allocate (varlist(max(d1max,ncand)),stat=st)
            if (st /= 0) go to 490
          end if
! The work array reals must be of size at least nfront.
          if (size(reals) < nfront) then
            deallocate (reals,stat=st)
            allocate (reals(int(nfront*multiplier)),stat=st)
            if (st /= 0) go to 490
          end if

          nfront_long = nfront
          lfa = (nfront_long*(nfront_long+nb+1))/2
          if (lfa > keep%lup) then
            info%flag = -30; return
          end if

          ilong = size_fa
          if (lfa > ilong) then
! write out contents of fa, allocate larger front and read back in
            locw = rtop
            lreals = cnvar
            ipp = ipc
            do j = 1,cnvar
              call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork,  &
                   locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,  &
                   keep%maxstore,keep%used)
!             write (6,*) 'ipp,ipp+lreals-1',ipp,ipp+lreals-1,j,cnvar
!             write (6,*) fa(ipp:ipp+lreals-1)
              if (flag < 0) go to 480
              ipp = ipp + lreals
              locw = locw + lreals
              lreals = lreals - 1
            end do
            rtop = rtop + length
            keep%rtopmx = max(keep%rtopmx,rtop)
            deallocate (fa,stat=st)
            lfa = ((nfront_long*multiplier)*((nfront_long*multiplier)+nb+1))/2
            if (lfa > keep%lup) lfa = (nfront_long*(nfront_long+nb+1))/2
            if (lfa > keep%lup) then
              info%flag = -30; return
            end if
            allocate (fa(lfa),stat=st)
            if (st /= 0) then
              info%flag = -30; info%stat = st; return
            end if
            size_fa = lfa
! read fa back in (no need to retain the temporary data)
            ipc = 1 + size_fa - length
            loc = rtop - length
            lreals = cnvar
            ipp = ipc
            do j = 1,cnvar
              call MA77_read_discard_real(keep%index(3),keep%rwork,loc,&
                   lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1)
!             write (6,*) 'ipp,ipp+lreals-1',ipp,ipp+lreals-1,j,cnvar
!             write (6,*) fa(ipp:ipp+lreals-1)
              if (flag < 0) go to 480
              ipp = ipp + lreals
              loc = loc + lreals
              lreals = lreals - 1
            end do
            rtop = rtop - length

          end if
!         write (6,*) 'node,nfront,nb,size_fa',node,nfront,nb,size_fa, &
!            real(nfront)*real(nfront)
        end if
! ip1 points to position after delayed columns
        ip1 = 1 + size_fa - (nvar_long*nvar_long+nvar_long)/2

! Read the variables for node from the main integer superfile.
        loc = keep%ifile(node)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
             rlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
        do i = 1,nvar
          ivar = abs(rlist(i))
          pos(ivar) = i
        end do

! Set ip(j) to point to start of col j in array fa.
        ip(1) = ip1
        do j = 2,nvar
          ip(j) = ip(j-1) + nvar - j + 2
        end do
!       write (6,*) 'ip',ip(1:nvar)

! read integer list for cnode from main integer file
        loc = keep%ifile(cnode) + keep%tree(cnode)%nelim
        call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
             varlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
! Set mapping from cnode into the front
        do j = 1,cnvar
          k = abs(varlist(j))
          map(j) = pos(k)
        end do

!       write (6,*) 'ic,splitp,non_leaf,cdelay',ic,splitp,non_leaf,cdelay
        if (ic == splitp .or. ic == non_leaf) then
          if (cnvar /= nvar) then
! Expand cnode into the front
            j1 = 1
            do j = 1,cnvar
              j2 = map(j)
! col. j of cnode maps to col. j2 of front
! set cols j1:j2-1 of front to zero
              fa(ip(j1):ip(j2)-1) = zero
! Expand col. j of cnode into col. j2 of front.
! col j of cnode is currently at col k=j+nvar-cnvar of front
              k = j + nvar - cnvar
              if (j2 /= k) then
                ipj2 = ip(j2)
                fa(ipj2:ip(j2+1)-1) = zero
                ipk = ip(k)
                do l = j,cnvar
                  i = map(l)
                  fa(ipj2+i-j2) = fa(ipk+l-j)
                end do
              else
! j2 = k, so all columns are now in place
                j1 = nvar + 1
                exit
              end if
              j1 = j2 + 1
            end do
! Check final cols are set to 0.
            klong = size_fa
            if (j1 <= nvar) fa(ip(j1):klong) = zero
          end if
        end if

        if (ic == splitp) then
! At split point.
! Merge children 1:splitp-1 into rows/columns 1:nfront of frontal matrix
          do 160 jc = splitp-1, 1, -1
            jcnode = keep%tree(node)%child(jc)
            nelim_jc = keep%tree(jcnode)%nelim
            loc = keep%ifile(jcnode) + nelim_jc
            cnvar = keep%size(jcnode) - nelim_jc
            call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
                 varlist,flag,keep%idata,-1)
            if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos (length n)
            do j = 1,cnvar
              k = abs(varlist(j))
              map(j) = pos(k)
            end do
! Read in reals for child from the top of stack, one col. at a time,
! using map to read directly into the front.
            lreals = cnvar
            length = cnvar
            length = (length*length + length)/2
            loc = rtop - length
! move pointer to top of stack
            rtop = loc
            do j = 1, cnvar
              ilong = map(j)
! start of column i+1 is
! k = ip1 + nvar*i - i*(i-1)/2
! so last entry in col. i is
! k - 1 = ip1 + nvar*i - i*(i-1)/2 - 1
! require length of part of we are reading into to be max(map) = nvar
! so start to read into k - 1 - nvar + 1 = k - nvar
! that is, ipp = k - nvar = ip1 + nvar*i - i*(i-1)/2 - nvar
              ipp = ip1 + nvar_long*(ilong-1) - (ilong*(ilong-1))/2
              call MA77_read_discard_real(keep%index(3), &
                   keep%rwork,loc,lreals,  &
                   fa(ipp:ipp+nvar-1),flag,keep%rdata,-1,map=map(j:cnvar))
              if (flag < 0) go to 480
              loc = loc + lreals
              lreals = lreals - 1
            end do
! End of loop over children
160       continue

! At this point fa holds the reals in the frontal matrix in packed form.
! If split point is not also last non-leaf child, write fa onto the stack
          if (ic /= non_leaf) then
            lreals = nvar
            loc = rtop
            ipp = ip1
            do j = 1,nvar
              call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork, &
                   loc,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1, &
                   keep%maxstore,keep%used)
              if (flag < 0) go to 480
              ipp = ipp + lreals
              loc = loc + lreals
              lreals = lreals - 1
            end do
            rtop = rtop + (nvar_long*nvar_long+nvar_long)/2
            keep%rtopmx = max(keep%rtopmx,rtop)
          end if

        else if (ic < non_leaf) then
! We have passed split point and not yet reached last non-leaf child.
! The reals for cnode are still in fa.
! Loop over the cols of cnode, merging the child into its parent.
            do j = 1,cnvar
              i = map(j)
! Read column i of the frontal matrix from the stack. Find its position.
               loc = rtop - (nvar_long*nvar_long+nvar_long)/2
! loc now points to the first entry in the front on the stack.
! Set loc to point to first entry in col. i and
! lreals to number of entries in col. i.
              ilong = i
              loc = loc + nvar_long*(ilong-1) - ((ilong-2)*(ilong-1))/2
              lreals = nvar - i + 1
              call MA77_read_real(keep%index(3),keep%rwork,loc,lreals,reals,&
                   flag,keep%rdata,-1)
              if (flag < 0) go to 480
! ipp points to start of col. j of cnode, which is held in fa
              jlong = j
              ipp = ipc + cnvar*(jlong-1) - ((jlong-2)*(jlong-1))/2
              do k = j,cnvar
                kk = map(k)
                reals(1+kk-i) = reals(1+kk-i) + fa(ipp+k-j)
              end do
! Write the updated reals back onto the stack
              call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork, &
                   loc,lreals,reals,flag,keep%rdata,-1,keep%maxstore,keep%used)
              if (flag < 0) go to 480
            end do

        else if (ic == non_leaf) then
! Dealing with last non-leaf child (and it is not split point)
! Read frontal matrix from top of stack, one column at a time and add
! in with last non-leaf child (which we have already expanded in fa).
          loc = rtop - (nvar_long*nvar_long+nvar_long)/2
          rtop = loc
          ipp = ip1
          lreals = nvar
          do j = 1,nvar
            call MA77_read_discard_real(keep%index(3),keep%rwork,loc,&
                 lreals,reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480
            fa(ipp:ipp+lreals-1) = fa(ipp:ipp+lreals-1) + reals(1:lreals)
            ipp = ipp + lreals
            loc = loc + lreals
            lreals = lreals - 1
          end do
!          lreals = (nvar*nvar+nvar)/2
!          call MA77_read_real(keep%index(3),keep%rwork,loc,lreals,&
!               fa(ip1:ip1+lreals-1),flag,keep%rdata,-1,add=.true.)

        end if

      if (info%flag >= 0) return

480   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

490   info%flag = -1
      info%stat = st

    end subroutine assemble64

!*************************

! End factorization

  end subroutine MA77_factor_solve_double

!*************************************************************************
   subroutine MA77_resid_double(nrhs,lx,x,lresid,resid,keep,control,info, &
      anorm_bnd)
! Compute the residual and optionally a bound on the norm of A.

    integer(short) :: nrhs
    integer(short) :: lx ! must be at least n
    integer(short) :: lresid ! must be at least n
    real(wp), intent(in) :: x(lx,nrhs)
    real(wp), intent(inout) :: resid(lresid,nrhs)
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info
    real(wp), intent(out),optional :: anorm_bnd

    real (wp), allocatable :: reals(:) !
    integer(short), allocatable :: varlist(:) ! used for variable lists
    integer(short), allocatable :: iwork(:) !
    real (wp), allocatable :: work(:) ! Work array (only if anorm_bnd present)

    real(wp) :: atemp
    integer(short) :: flag
    integer(short) :: i ! temporary variable
    integer(short) :: ielt
    integer(short) :: inelrs
    integer(short) :: j
    integer(short) :: k
    integer(short) :: l
    integer(short) :: lfa
    integer(short) :: ll
    integer(short) :: mvar
    integer(short) :: n
    integer(short) :: nout
    integer(short) :: nvar
    integer(short) :: st ! stat parameter
    integer(long) :: loc,locr

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -24   Error in size of x
! -25   Error in size of resid

! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
     ' Entering MA77_resid with:'
      write (control%unit_diagnostics,'(a,4(/a,i12),(/a,i12/a,i12/a,i12))') &
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
       nrhs, &
     ' lx                                                      = ', &
       lx, &
     ' lresid                                                  = ', &
       lresid
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_resid'

! Check status parameter
    if (keep%status /= 3) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    info%flag = keep%flag
    n = keep%n
    if (n == 0) then
      if (present(anorm_bnd)) anorm_bnd = zero
      return
    end if

    if (lx < n) then
      info%flag = -24
      keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag)
      if (nout >= 0) write (nout,'(a,i8,a,i8)') &
    ' Increase lx from ', lx, ' to at least ', n
      return
     end if

    if (nrhs < 1) then
      info%flag = -24
      keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag)
      if (nout >= 0) write (nout,'(a,i8,a,i8)') &
    ' nrhs must be at least 1. nrhs = ', nrhs
      return
     end if

    if (lresid < n) then
      info%flag = -25
      call MA77_print_iflag(keep,nout,info%flag)
      if (nout >= 0) write (nout,'(a,i8,a,i8)') &
   ' Increase lresid from ', lresid, ' to at least ', n
      return
    end if

      mvar = keep%mvar
      if (keep%element_input) then
        deallocate (varlist,stat=st)
        allocate (varlist(mvar),stat=st)
        if (st /= 0) go to 95
        if (keep%index(2) >= 0) then
          lfa = (mvar*mvar + mvar)/2
          deallocate (reals,stat=st)
          allocate (reals(lfa),stat=st)
          if (st /= 0) go to 95
        end if

        if (present(anorm_bnd)) then
! norm of A is required so need extra array
          deallocate (work,stat=st)
          allocate (work(n),stat=st)
          if (st /= 0) go to 95
          work = zero
          do ielt = 1,keep%nelt
            nvar = abs(keep%size(ielt))
            if (nvar == 0) cycle
            loc = keep%ifile(ielt)
            if (keep%index(1) >= 0) then
! read in variable list
              call of01_read(keep%index(1),loc,nvar,varlist,flag,keep%idata, &
                   lp=-1)
              if (flag < 0) go to 97
            else
! Direct access files not in use.
              varlist(1:nvar) = keep%imain(loc:loc+nvar-1)
            end if
! read in the reals
            inelrs = (nvar*nvar+nvar)/2
            locr = keep%rfile(ielt) - inelrs
            if (keep%index(2) >= 0) then
              call of01_read(keep%index(2),locr,inelrs,reals,flag, &
                   keep%rdata,lp=-1)
              if (flag < 0) go to 97
              k = 1
              do i = 1,nvar
                l = varlist(i)
                atemp = -reals(k)
                resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(l,1:nrhs)
                work(l) = work(l) + abs(atemp)
                k = k + 1
                do j = i+1,nvar
                  ll = varlist(j)
                  atemp = -reals(k)
                  resid(ll,1:nrhs) = resid(ll,1:nrhs) + atemp*x(l,1:nrhs)
                  resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(ll,1:nrhs)
                  atemp = abs(atemp)
                  work(l) = work(l) + atemp
                  work(ll) = work(ll) + atemp
                  k = k + 1
                end do
              end do
            else
! Direct access files not in use
              k = 1
              do i = 1,nvar
                l = varlist(i)
                atemp = -keep%rmain(locr+k-1)
                resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(l,1:nrhs)
                work(l) = work(l) + abs(atemp)
                k = k + 1
                do j = i+1,nvar
                  ll = varlist(j)
                  atemp = -keep%rmain(locr+k-1)
                  resid(ll,1:nrhs) = resid(ll,1:nrhs) + atemp*x(l,1:nrhs)
                  resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(ll,1:nrhs)
                  atemp = abs(atemp)
                  work(l) = work(l) + atemp
                  work(ll) = work(ll) + atemp
                  k = k + 1
                end do
              end do
            end if
          end do
          anorm_bnd = work(1)
          do i = 2,n
            anorm_bnd = max(anorm_bnd,work(i))
          end do
        else
! norm of A not required
          do ielt = 1,keep%nelt
            nvar = abs(keep%size(ielt))
            if (nvar == 0) cycle
            loc = keep%ifile(ielt)
! read in variable list
            if (keep%index(1) >= 0) then
              call of01_read(keep%index(1),loc,nvar,varlist,flag, &
                   keep%idata,lp=-1)
              if (flag < 0) go to 97
            else
              varlist(1:nvar) = keep%imain(loc:loc+nvar-1)
            end if
! read in the reals
            inelrs = (nvar*nvar+nvar)/2
            locr = keep%rfile(ielt) - inelrs
            if (keep%index(2) >= 0) then
              call of01_read(keep%index(2),locr,inelrs,reals,flag, &
                   keep%rdata,lp=-1)
              if (flag < 0) go to 97
              k = 1
              do i = 1,nvar
                l = varlist(i)
                atemp = -reals(k)
                resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(l,1:nrhs)
                k = k + 1
                do j = i+1,nvar
                  ll = varlist(j)
                  atemp = -reals(k)
                  resid(ll,1:nrhs) = resid(ll,1:nrhs) + atemp*x(l,1:nrhs)
                  resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(ll,1:nrhs)
                  k = k + 1
                end do
              end do
            else
              k = 1
              do i = 1,nvar
                l = varlist(i)
                atemp = -keep%rmain(locr+k-1)
                resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(l,1:nrhs)
                k = k + 1
                do j = i+1,nvar
                  ll = varlist(j)
                  atemp = -keep%rmain(locr+k-1)
                  resid(ll,1:nrhs) = resid(ll,1:nrhs) + atemp*x(l,1:nrhs)
                  resid(l,1:nrhs) = resid(l,1:nrhs) + atemp*x(ll,1:nrhs)
                  k = k + 1
                end do
              end do
            end if
          end do
        end if
      else
! row entry
        if (keep%index(1) >= 0) then
          deallocate (varlist,stat=st)
          allocate (varlist(mvar),stat=st)
          if (st /= 0) go to 95
        end if
        if (keep%index(2) >= 0) then
          deallocate (reals,stat=st)
          allocate (reals(mvar),stat=st)
          if (st /= 0) go to 95
        end if
        if (present(anorm_bnd)) then
! norm of A is required so need extra array
          deallocate (work,stat=st)
          allocate (work(n),stat=st)
          if (st /= 0) go to 95
          work = zero
          do ielt = 1,keep%nelt
            nvar = abs(keep%size(ielt))
            if (nvar == 0) cycle
            loc = keep%ifile(ielt)
            locr = keep%rfile(ielt) - nvar
! read in variable list
            if (keep%index(1) >= 0) then
              call of01_read(keep%index(1),loc,nvar,varlist,flag, &
                   keep%idata,lp=-1)
              if (flag < 0) go to 97
! read in the reals
              if (keep%index(2) >= 0) then
                call of01_read(keep%index(2),locr,nvar,reals,flag, &
                     keep%rdata,lp=-1)
                if (flag < 0) go to 97
                do i = 1,nvar
                  l = varlist(i)
                  atemp = -reals(i)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                  work(l) = work(l) + abs(atemp)
                end do
              else
                do i = 1,nvar
                  l = varlist(i)
                  atemp = -keep%rmain(locr+i-1)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                  work(l) = work(l) + abs(atemp)
                end do
              end if
            else
              if (keep%index(2) >= 0) then
                call of01_read(keep%index(2),locr,nvar,reals,flag, &
                     keep%rdata,lp=-1)
                if (flag < 0) go to 97
                do i = 1,nvar
                  l = keep%imain(loc+i-1)
                  atemp = -reals(i)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                  work(l) = work(l) + abs(atemp)
                end do
              else
                do i = 1,nvar
                  l = keep%imain(loc+i-1)
                  atemp = -keep%rmain(locr+i-1)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                  work(l) = work(l) + abs(atemp)
                end do
              end if
            end if
          end do
          anorm_bnd = work(1)
          do i = 2,n
            anorm_bnd = max(anorm_bnd,work(i))
          end do
        else
! norm of A not required
          do ielt = 1,keep%nelt
            nvar = abs(keep%size(ielt))
            if (nvar == 0) cycle
            loc = keep%ifile(ielt)
            locr = keep%rfile(ielt) - nvar
! read in variable list
            if (keep%index(1) >= 0) then
              call of01_read(keep%index(1),loc,nvar,varlist,flag,keep%idata, &
                             lp=-1)
              if (flag < 0) go to 97
! read in the reals
              if (keep%index(2) >= 0) then
                call of01_read(keep%index(2),locr,nvar,reals,flag, &
                     keep%rdata,lp=-1)
                if (flag < 0) go to 97
                do i = 1,nvar
                  l = varlist(i)
                  atemp = -reals(i)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                end do
              else
                do i = 1,nvar
                  l = varlist(i)
                  atemp = -keep%rmain(locr+i-1)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                end do
              end if
            else
! integer superfile not in use
              if (keep%index(2) >= 0) then
                call of01_read(keep%index(2),locr,nvar,reals,flag, &
                     keep%rdata,lp=-1)
                if (flag < 0) go to 97
                do i = 1,nvar
                  l = keep%imain(loc+i-1)
                  atemp = -reals(i)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                end do
              else
                do i = 1,nvar
                  l = keep%imain(loc+i-1)
                  atemp = -keep%rmain(locr+i-1)
                  resid(ielt,1:nrhs) = resid(ielt,1:nrhs) + atemp*x(l,1:nrhs)
                end do
              end if
            end if
          end do
        end if

      end if
    if (present(anorm_bnd)) then
      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
        write (control%unit_diagnostics,'(/a,es12.4)') &
      ' Leaving MA77_resid with        anorm_bnd                = ', anorm_bnd
    end if

    go to 99

 95 info%flag = -1
    info%stat = st
    call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
    go to 99

 97 info%iostat = keep%idata%iostat
    info%stat = keep%idata%stat
    info%flag = flag
    if (flag == -17) info%flag = -16
    call MA77_print_iflag(keep,nout,info%flag,st=info%stat, &
         ios=keep%idata%iostat)

 99 deallocate (varlist,stat=st)
    deallocate (iwork,stat=st)
    deallocate (work,stat=st)
    deallocate (reals,stat=st)

    end subroutine MA77_resid_double
!*************************************************************************

  subroutine MA77_solve_double(nrhs,lx,x,keep,control,info,scale,job)

! Solve phase. Optionally performs only the forward or backward sub.
!    use indef
    integer(short), intent (in) :: nrhs
    integer(short), intent (in) :: lx
    real (wp), intent (inout) :: x(lx,nrhs) ! On entry, x must
!             be set so that if i has been used to index a variable,
!             x(i,j) is the corresponding component of the
!             right-hand side for the jth system (j = 1,2,..., nrhs).
!             On exit, if i has been used to index a variable,
!             x(i,j) holds solution for variable i to system j
! For details of keep, control, info : see derived type description
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info
    real(wp), intent(in), optional :: scale(:) ! if present,
!      must hold scaling factors. If scale was present on
!      call to ma77_factor , then scale should also be present here.
   integer(short), optional, intent (in) :: job  ! used to indicate whether
!             partial solution required
!             job = 1 : forward eliminations only (PLX = B)
!             job = 2 : diagonal solve (DX = B) (indefinite case only)
!             job = 3 : backsubs only ((PL)^TX = B)
!             job = 4 : diag and backsubs (D(PL)^TX = B) (indefinite case only)
!             job absent: complete solve performed

    real (wp), allocatable :: xlocal(:,:) ! temporary dense vector for
!              performing forward eliminations and back substitutions.
    real (wp), allocatable :: d(:)   ! used to hold diagonal entries
    real (wp), allocatable :: fa(:)  ! used to hold factor entries
    real (wp), allocatable :: w54(:) ! workspace for ma54 (nrhs > 4)
    integer(short), allocatable :: varlist(:) ! used for variable lists
    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -20   job is out of range
! -24   Error in size of x
! -39   size(scale) too small or scale not present when expected or
!       present when not expected

    integer(short) :: depth ! depth in tree
    integer(short) :: dreals ! (indefinite case only)
    integer(short) :: flag ! local error flag
    integer(short) :: ir ! do loop variable
    integer(short) :: i ! temporary variable
    integer(short) :: j ! temporary variable
    integer(short) :: jvar ! variable belonging to node
    integer(long)  :: lq ! number of factor entries to read (indefinite)
    integer(long)  :: loc ! location in file
    integer(short) :: local_job ! local job parameter
    integer(short) :: lreals ! number of reals to be read at node
    integer(short) :: lw54 ! length of workarray w54
    integer(long)  :: lcopy ! temporary copy of locend
    integer(long)  :: locend ! initialised to keep%rfree (first free position
!          after end of factor storage)
    integer(long)  :: locendi ! initialised to keep%ifree (first free position
!          after end of factor storage) (integer data)
    integer(long)  :: locfac ! initialised to start of factor storage
    integer(long)  :: locint ! initialised to start of integer factor
    integer(short) :: maxlen ! holds largest number of variables in
!          an element list
    integer(short) :: nchild ! number of children of node in tree
    integer(short) :: nelim ! number of eliminations at node
    integer(long)  :: nelim_long ! number of eliminations at node
    integer(short) :: nfront ! number of variables in node
    integer(long)  :: nfront_long ! number of variables in node
    integer(short) :: nb ! set to keep%nb
    integer(short) :: n ! set to keep%n
    integer(short) :: nelt ! set to keep%nelt
    integer(short) :: node ! node in tree
    integer(short) :: nout ! output unit
    integer(short) :: nroot ! number of roots (components)
    integer(short) :: pnode ! parent node in tree
    integer(long)  :: posdiag ! initialised to start of diagonal entries
!          storage. Incremented as integers read (indefinite case only)
    integer(short) :: root ! root of tree
    integer(short) :: st ! stat parameter
    logical :: ldiag_solve ! controls whether a diagonal solve
!          is performed on the call to the backsubstitution.
!          it is used by the subroutine back in the indefinite case.
!          If set to false, then diag. solve not performed (job = 3).
!          If set to true, then diag. solve is performed.

!    tblas = 0.0; sizes = 0

! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
     ' Entering MA77_solve with:'
      write (control%unit_diagnostics,'(a,4(/a,i12),(/a,i12/a,i12))') &
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
       nrhs, &
     ' lx                                                      = ', &
       lx
    end if
    if (present(scale)) then
      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
      write (control%unit_diagnostics,'(a)') &
     ' Scaling factors supplied. '
    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_solve'

    info%flag = keep%flag
    info%nio_read(1)  = keep%idata%nio_read
    info%nio_read(2)  = keep%rdata%nio_read
    info%nio_write(1) = keep%idata%nio_write
    info%nio_write(2) = keep%rdata%nio_write
    info%nwd_read(1)  = keep%idata%nwd_read
    info%nwd_read(2)  = keep%rdata%nwd_read
    info%nwd_write(1) = keep%idata%nwd_write
    info%nwd_write(2) = keep%rdata%nwd_write

! Check status parameter
    if (keep%status /= 3) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    n = keep%n
    if (keep%n == 0) return
    nelt = keep%nelt

    if (lx < n) then
      info%flag = -24
      keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag)
      if (nout >= 0) write (nout,'(a,i8,a,i8)') &
    ' Increase lx from ', lx, ' to at least ', n
      return
     end if

    if (nrhs < 1) then
      info%flag = -24
      keep%status = -1
      call MA77_print_iflag(keep,nout,info%flag)
      if (nout >= 0) write (nout,'(a,i8,a,i8)') &
    ' nrhs must be at least 1. nrhs = ', nrhs
      return
    end if

    if (present(scale)) then
      if (keep%scale == 0) info%flag = -39      
      if (size(scale) < n) info%flag = -39
    else
      if (keep%scale == 1) info%flag = -39
    end if
    if (info%flag == -39) then
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    local_job = 0
    if (present(job)) then
      if (job < 1 .or. job > 4) info%flag = -20
      if (keep%pos_def .and. job == 2) info%flag = -20
      if (keep%pos_def .and. job == 4) info%flag = -20
      if (info%flag == -20) then
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if
      local_job = job
    end if

    deallocate (varlist,stat=st)
    deallocate (xlocal,stat=st)
    deallocate (w54,stat=st)
    deallocate (d,stat=st)
    deallocate (fa,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)

    allocate (count(0:keep%maxdepth),cnode(0:keep%maxdepth),stat=st)
    if (st /= 0) goto 95

    maxlen = keep%maxlen
    if (.not. keep%pos_def) maxlen = max(maxlen,keep%maxfrontb)

    nb = keep%nb
    if (keep%pos_def) then
      lw54 = 1
      if (nrhs >=4) lw54 = 2*nb*maxlen
      allocate (w54(lw54),stat=st)
      if (st /= 0) go to 95
    end if

    if (keep%index(2) >= 0) then
      if (local_job /= 2) then
        allocate (fa(keep%maxfa),stat=st)
        if (st /= 0) go to 95
      end if
      if (.not. keep%pos_def .and. local_job /= 1) then
        allocate (d(2*keep%maxelim_actual),stat=st)
        if (st /= 0) go to 95
      end if
    end if

    if (keep%index(1) >= 0) then
      allocate (varlist(maxlen+1),stat=st)
      if (st /= 0) go to 95
    end if
    allocate (xlocal(maxlen,nrhs),stat=st)
    if (st /= 0) go to 95

    nroot = size(keep%roots)

    if (present(scale)) then
      if (local_job == 0 .or. local_job == 1) then
        do i = 1,n
          x(i,1:nrhs) = scale(i)*x(i,1:nrhs)
        end do
      end if
    end if

    ldiag_solve = .true.

    if (local_job == 0) then
! full system solve.
! Forward eliminations. Initialise locfac/locint to the start of
! the factor storage (skip diagonal entries that are stored first)
      locfac = keep%posfac + n
      if (.not. keep%pos_def) locfac = locfac + n
      locint = keep%posint

    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. 
        call forward(node)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 96
    end do ! roots


! At this point, array x holds partial solution vector(s).
! It will be overwritten by solution.
! Back substitutions. Initialise locend/locendi to first free position
! after the end of the factor storage (we must read
! records in reverse order to the written order)
      locend  = keep%rfree
      locendi = keep%ifree
      posdiag = keep%dfree

      do ir = nroot, 1, -1
        root = keep%roots(ir)
! Visit each node in tree rooted at root in reverse depth-first search order
        node = root
        cnode(0) = root; count(0) = 0
        depth = 1
        count(1) = 0 ! number of visited children
        do
          cnode(depth) = node
          if (count(depth) == 0 .and. node > nelt) then
            call back(node)
            if (info%flag < 0) exit
            nchild = size(keep%tree(node)%child)
            if (nchild == 0) then
              if (node == root) exit
            else
! Descend to last child of node 
              count(depth) = 1
              node = keep%tree(node)%child(nchild)
              depth = depth + 1
              count(depth) = 0
              cycle
            end if
          end if

! Parent node
          pnode = cnode(depth-1)
          nchild = size(keep%tree(pnode)%child)
          if (count(depth-1) < nchild) then
! There is another sibling of node  - go to it
            count(depth-1) = count(depth-1) + 1
            node = keep%tree(pnode)%child(nchild-count(depth-1)+1)
            count(depth)=0
            cycle
          end if

! No siblings remain so go to parent
          depth = depth - 1
          node = pnode
! If we have reached the root, we are done
          if (depth <= 1) exit

        end do
        if (info%flag < 0) go to 96
      end do

    else if (local_job == 1) then
! Foward sub.
      locfac = keep%posfac + n
      if (.not. keep%pos_def) locfac = locfac + n
      locint = keep%posint
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. 
        call forward(node)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 96
    end do ! roots

    else if (local_job == 2) then
! Diag. solve
      locint = keep%posint
      posdiag = keep%dfree
      posdiag = keep%posfac
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node.
        call diag_solve(node)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 96
    end do ! roots

    else if (local_job == 3 .or. local_job == 4) then
! Back sub. (without diagonal solve in indefinite case if job = 3)
      locend  = keep%rfree
      locendi = keep%ifree
      posdiag = keep%dfree
      if (local_job == 3) ldiag_solve = .false.
      do ir = nroot, 1, -1
        root = keep%roots(ir)
! Visit each node in tree rooted at root in reverse depth-first search order
        node = root
        cnode(0) = root; count(0) = 0
        depth = 1
        count(1) = 0 ! number of visited children
        do
          cnode(depth) = node
          if (count(depth) == 0 .and. node > nelt) then
            call back(node)
            if (info%flag < 0) exit
            nchild = size(keep%tree(node)%child)
            if (nchild == 0) then
              if (node == root) exit
            else
! Descend to last child of node 
              count(depth) = 1
              node = keep%tree(node)%child(nchild)
              depth = depth + 1
              count(depth) = 0
              cycle
            end if
          end if

! Parent node
          pnode = cnode(depth-1)
          nchild = size(keep%tree(pnode)%child)
          if (count(depth-1) < nchild) then
! There is another sibling of node  - go to it
            count(depth-1) = count(depth-1) + 1
            node = keep%tree(pnode)%child(nchild-count(depth-1)+1)
            count(depth)=0
            cycle
          end if

! No siblings remain so go to parent
          depth = depth - 1
          node = pnode
! If we have reached the root, we are done
          if (depth <= 1) exit

        end do
        if (info%flag < 0) go to 96
      end do

    end if

    if (present(scale)) then
      if (local_job == 0 .or. local_job == 3 .or. local_job == 4) then
        do i = 1,n
          x(i,1:nrhs) = scale(i)*x(i,1:nrhs)
        end do
      end if
    end if

    info%nio_read(1)  = keep%idata%nio_read  - info%nio_read(1)
    info%nio_read(2)  = keep%rdata%nio_read  - info%nio_read(2)
    info%nio_write(1) = keep%idata%nio_write - info%nio_write(1)
    info%nio_write(2) = keep%rdata%nio_write - info%nio_write(2)
    info%nwd_read(1)  = keep%idata%nwd_read  - info%nwd_read(1)
    info%nwd_read(2)  = keep%rdata%nwd_read  - info%nwd_read(2)
    info%nwd_write(1) = keep%idata%nwd_write - info%nwd_write(1)
    info%nwd_write(2) = keep%rdata%nwd_write - info%nwd_write(2)

    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(/a)') &
     ' Completed solve with:'
      write (control%unit_diagnostics, &
     '(a,(/a,i12),4(/a/2es12.4))') &
     ' information parameters (info%) :', &
     ' flag                Error flag                          = ', &
       info%flag, &
  ' nio_read(1:2)       Number of records read from disk by OF01_read   = ',&
       real(info%nio_read(1:2)), &
  ' nio_write(1:2)      Number of records written to disk by OF01_write = ',&
       real(info%nio_write(1:2)), &
     ' nwd_read(1:2)       Number of scalars read by OF01_read     = ', &
       real(info%nwd_read(1:2)), &
     ' nwd_write(1:2)      Number of scalars written by OF01_write = ', &
       real(info%nwd_write(1:2))
    end if

!    write (6,'(a,f10.3)') 'tblas = ',tblas
!    write (6,'(a,10i7)') 'sizes =',sizes(1:10)

    if (info%flag >= 0) go to 100

 95   info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      go to 100

 96   call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)

! Deallocate arrays we have done with.
 100  deallocate (varlist,stat=st)
      deallocate (xlocal,stat=st)
      deallocate (fa,stat=st)
      deallocate (d,stat=st)
      deallocate (w54,stat=st)
      deallocate (count,stat=st)
      deallocate (cnode,stat=st)
  contains

!*******************************

    subroutine forward(node)
! This subroutine performs forward eliminations at node.
! User-supplied rhs vectors are overwritten by the partial solution.

      integer(short), intent (in) :: node

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

      if (keep%pos_def) then
        nelim = keep%tree(node)%nelim
        nfront = abs(keep%size(node))
        loc = keep%ifile(node)
        if (keep%index(1) >= 0) then
! Read integer data from file
          call of01_read(keep%index(1),loc,nfront,varlist,flag, &
               keep%idata,lp=-1)
          if (flag < 0) go to 480
! Copy part of right-hand side array x corresponding to the variables
! that are in varlist into local array xlocal
          do i = 1, nfront
            jvar = varlist(i)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        else
! Integers are held in keep%imain
          do i = 1, nfront
            jvar = keep%imain(loc+i-1)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        end if
!       write (6,'(6es12.4)') xlocal(1:nfront,1)

        if (keep%index(2) >= 0) then
! Read in required part of factor and do forward elimination ops.
          lreals = nfront
          i = 1
          do j = 1,nelim
            call of01_read(keep%index(2),locfac,lreals,fa(i:i+lreals-1), &
                 flag,keep%rdata,lp=-1)
            if (flag < 0) go to 485
            locfac = locfac + lreals
            i = i + lreals
            lreals = lreals - 1
          end do
          if (nrhs >= 4) then
            call ma54_forward2(nfront,nelim,nb,nrhs,fa,xlocal,maxlen,nb, &
                 w54,flag)
          else
            do j = 1,nrhs
              call ma54_forward1(nfront,nelim,nb,fa,xlocal(:,j),flag)
            end do
          end if
        else
! Factor entries are held in keep%rmain
          nfront_long = nfront
          nelim_long = nelim
          lq = (nelim_long*nelim_long+nelim_long)/2 + &
                    nelim_long*(nfront_long-nelim_long)
          if (nrhs >= 4) then
            call ma54_forward2(nfront,nelim,nb,nrhs, &
                 keep%rmain(locfac:locfac+lq-1),xlocal,maxlen,nb,w54,flag)
          else
            do j = 1,nrhs
              call ma54_forward1(nfront,nelim,nb, &
                   keep%rmain(locfac:locfac+lq-1),xlocal(:,j),flag)
            end do
          end if
          locfac = locfac + lq
        end if
!        if (flag /= 0) then
!          write (6,*) 'error forward: flag=', flag
!        end if

!       write (6,'(a,6es12.4)') 'after  ',xlocal(1:nelim,1)

      else
! Indefinite case. The variable lists are held after the main integer lists
        nfront = keep%size_ind(node)
        if (keep%index(1) >= 0) then
          call of01_read(keep%index(1),locint,nfront+1,varlist,flag, &
               keep%idata,lp=-1)
          if (flag < 0) go to 480
! Increment locint so that it points to start of integer data for next node
          locint = locint + nfront + 1
          nelim = varlist(nfront+1)
          if (nelim == 0) return
! Copy part of right-hand side array x corresponding to the variables
! that are in varlist into local array xlocal
          do i = 1, nfront
            jvar = varlist(i)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        else
! File not in use. Integers are held in keep%imain
          loc = locint
          locint = locint + nfront + 1
          nelim = keep%imain(loc+nfront)
          if (nelim == 0) return
          do i = 1, nfront
            jvar = keep%imain(loc+i-1)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        end if

        lq = nelim
        lq = (lq*lq+lq)/2 + lq*(nfront-lq)
        if (keep%index(2) >= 0) then
! read required part of factor and do forward eliminations
          lreals = nfront
          i = 1
          do j = 1,nelim
            call of01_read(keep%index(2),locfac,lreals,fa(i:i+lreals-1), &
                 flag,keep%rdata,lp=-1)
            if (flag < 0) go to 485
            locfac = locfac + lreals
            i = i + lreals
            lreals = lreals - 1
          end do

          if (nrhs >= 4) then
            call ma64_solveL2(nfront,nelim,nb,nrhs,xlocal,maxlen,flag,fa,lq)
          else
            do j = 1,nrhs
              call ma64_solveL1(nfront,nelim,nb,xlocal(:,j),flag,fa,lq)
             end do
          end if
        else
! factors are in keep%rmain
          if (nrhs >= 4) then
            call ma64_solveL2(nfront,nelim,nb,nrhs,xlocal,maxlen, &
                 flag,keep%rmain(locfac:locfac+lq-1),lq)
          else
            do j = 1,nrhs
              call ma64_solveL1(nfront,nelim,nb,xlocal(:,j),flag, &
                   keep%rmain(locfac:locfac+lq-1),lq)
            end do
          end if
          locfac = locfac + lq
        end if
!       write (6,'(a,6es12.4)') 'after  ',xlocal(1:nfront,1)
      end if

! copy xlocal back to x (we overwrite the given right-hand
! side array with the partial solution)
      if (keep%index(1) >= 0) then
        do i = 1, nfront
          jvar = varlist(i)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
        end do
      else
        do i = 1, nfront
          jvar = keep%imain(loc+i-1)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
        end do
      end if
!     write (6,'(a,6es12.4)') 'after  ',x(1:keep%n,1)

      if (info%flag >= 0) go to 500

480   info%iostat = keep%idata%iostat
      info%stat = keep%idata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      go to 500

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      go to 500

500   return

    end subroutine forward

!*******************************
    subroutine back(node)
! This subroutine performs backsubstitution on
! the subtree rooted at node.

      integer(short), intent (in) :: node
! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

      if (keep%pos_def) then
        nelim = keep%tree(node)%nelim
        nfront = keep%size(node)
        loc = keep%ifile(node)
        if (keep%index(1) >= 0) then
          call of01_read(keep%index(1),loc,nfront,varlist,flag,&
               keep%idata,lp=-1)
          if (flag < 0) go to 480
! copy partial solution which is held in x into local array
          do i = 1, nfront
            jvar = varlist(i)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        else
          do i = 1, nfront
            jvar = keep%imain(loc+i-1)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        end if

        nfront_long = nfront
        nelim_long = nelim
        lq = (nelim_long*nelim_long+nelim_long)/2 + &
                  nelim_long*(nfront_long-nelim_long)
        locend = locend - lq
        if (keep%index(2) >= 0) then
! read required part of factor and do back sub. ops.
          lreals = nfront
          i = 1
          do j = 1,nelim
            call of01_read(keep%index(2),locend,lreals,fa(i:i+lreals-1), &
                flag,keep%rdata,lp=-1)
            if (flag < 0) go to 485
            locend = locend + lreals
            i = i + lreals
            lreals = lreals - 1
          end do
          locend = locend - lq
!         write (6,'(6es12.4)') xlocal(1:nelim,1)
          if (nrhs >= 4) then
            call ma54_back2(nfront,nelim,nb,nrhs,fa,xlocal,maxlen,nb,&
                 w54,flag)
          else
            do j = 1,nrhs
              call ma54_back1(nfront,nelim,nb,fa,xlocal(:,j),flag)
            end do
          end if
        else
          if (nrhs >= 4) then
            call ma54_back2(nfront,nelim,nb,nrhs, &
                 keep%rmain(locend:locend+lq-1),xlocal,maxlen,nb,w54,flag)
          else
            do j = 1,nrhs
              call ma54_back1(nfront,nelim,nb, &
                   keep%rmain(locend:locend+lq-1),xlocal(:,j),flag)
            end do
          end if
       end if

      else

! In indefinite case, lists are held after main initial integer lists.
        nfront = keep%size_ind(node)
        locendi = locendi - nfront - 1
        if (keep%index(1) >= 0) then
          call of01_read(keep%index(1),locendi,nfront+1,varlist,flag,&
               keep%idata,lp=-1)
          if (flag < 0) go to 480
          nelim = varlist(nfront+1)
          if (nelim == 0) go to 100
! copy partial solution which is held in x into local array
!         write (6,*) 'partial solution'
          do i = 1, nfront
            jvar = varlist(i)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        else
! file not in use
          nelim = keep%imain(locendi+nfront)
          if (nelim == 0) go to 100
          loc = locendi
! copy partial solution which is held in x into local array
!         write (6,*) 'partial solution'
          do i = 1, nfront
            jvar = keep%imain(loc+i-1)
            xlocal(i,1:nrhs) = x(jvar,1:nrhs)
          end do
        end if

        lq = nelim
        lq = (lq*lq+lq)/2 + lq*(nfront-lq)
        locend = locend - lq

        if (ldiag_solve) then
! Diagonal solve is to be performed

          dreals = 2*nelim
          posdiag = posdiag - dreals
          if (keep%index(2) >= 0) then
            lcopy = locend
            lreals = nfront
            i = 1
            do j = 1,nelim
              call of01_read(keep%index(2),locend,lreals,fa(i:i+lreals-1), &
                   flag,keep%rdata,lp=-1)
              if (flag < 0) go to 485
              locend = locend + lreals
              i = i + lreals
              lreals = lreals - 1
            end do
! reset locend
            locend = lcopy

! read in diagonal entries
            call of01_read(keep%index(2),posdiag,dreals, &
                 d,flag,keep%rdata,lp=-1)
            if (flag < 0) go to 485

            if (nrhs >= 4) then
              call ma64_solveDLT2(nfront,nelim,nb,nrhs,xlocal,maxlen, &
                   flag,fa,lq,d)
            else
              do j = 1,nrhs
                call ma64_solveDLT1(nfront,nelim,nb,xlocal(:,j),flag,fa,lq,d)
              end do
            end if

          else
            if (nrhs >= 4) then
              call ma64_solveDLT2(nfront,nelim,nb,nrhs,xlocal,maxlen, &
                   flag,keep%rmain(locend:locend+lq-1),lq, &
                   keep%rmain(posdiag:posdiag+dreals-1))
            else
              do j = 1,nrhs
                call ma64_solveDLT1(nfront,nelim,nb,xlocal(:,j),flag, &
                   keep%rmain(locend:locend+lq-1),lq, &
                   keep%rmain(posdiag:posdiag+dreals-1))
              end do
            end if
          end if

        else
! diagonal solve not wanted (this is on a job=3 call to ma77_solve)
          if (keep%index(2) >= 0) then
            lcopy = locend
            lreals = nfront
            i = 1
            do j = 1,nelim
              call of01_read(keep%index(2),locend,lreals,fa(i:i+lreals-1), &
                   flag,keep%rdata,lp=-1)
              if (flag < 0) go to 485
              locend = locend + lreals
              i = i + lreals
              lreals = lreals - 1
            end do
! reset locend
            locend = lcopy

            if (nrhs >= 4) then
              call ma64_solveLT2(nfront,nelim,nb,nrhs,xlocal,maxlen, &
                   flag,fa,lq)
            else
              do j = 1,nrhs
                call ma64_solveLT1(nfront,nelim,nb,xlocal(:,j),flag,fa,lq)
              end do
            end if

          else
            if (nrhs >= 4) then
              call ma64_solveLT2(nfront,nelim,nb,nrhs,xlocal,maxlen, &
                   flag,keep%rmain(locend:locend+lq-1),lq)
            else
              do j = 1,nrhs
                call ma64_solveLT1(nfront,nelim,nb,xlocal(:,j),flag, &
                   keep%rmain(locend:locend+lq-1),lq)
              end do
            end if
          end if

        end if

      end if

! copy part of xlocal corresponding to eliminated variables back into x
      if (keep%index(1) >= 0) then
        do i = 1, nelim
          jvar = varlist(i)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
!         write (6,*) jvar,x(jvar,1)
        end do
      else
        do i = 1, nelim
          jvar = keep%imain(loc+i-1)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
!         write (6,*) jvar,x(jvar,1)
        end do
      end if
!     write (6,'(6es12.4)') xlocal(1:nelim,1)


 100  if (info%flag >= 0) go to 500

480   info%iostat = keep%idata%iostat
      info%stat = keep%idata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      go to 500

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      go to 500

500   return

    end subroutine back
!*******************************

    subroutine diag_solve(node)
! This subroutine performs diagonal solve at node.

      integer(short), intent (in) :: node

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!      for all path)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

! The variable lists are held after the main integer lists
      nfront = keep%size_ind(node)
      if (keep%index(1) >= 0) then
        call of01_read(keep%index(1),locint,nfront+1,varlist,flag, &
             keep%idata,lp=-1)
        if (flag < 0) go to 480
! Increment locint so that it points to start of integer data for next node
        locint = locint + nfront + 1
        nelim = varlist(nfront+1)
        if (nelim == 0) return
! Copy part of right-hand side array x corresponding to the variables
! that are in varlist into local array xlocal
        do i = 1, nfront
          jvar = varlist(i)
          xlocal(i,1:nrhs) = x(jvar,1:nrhs)
        end do
      else
! File not in use. Integers are held in keep%imain
        loc = locint
        locint = locint + nfront + 1
        nelim = keep%imain(loc+nfront)
        if (nelim == 0) return
!! Copy into varlist
!!      varlist(1:nelim) = keep%imain(loc:loc+nelim-1)
        do i = 1, nfront
          jvar = keep%imain(loc+i-1)
          xlocal(i,1:nrhs) = x(jvar,1:nrhs)
        end do
      end if

      dreals = 2*nelim
      if (keep%index(2) >= 0) then
! read required part of d
        call of01_read(keep%index(2),posdiag,dreals,d, &
             flag,keep%rdata,lp=-1)
        if (flag < 0) go to 485
        if (nrhs >= 4) then
          call ma64_solveD2(nfront,nelim,nrhs,xlocal,maxlen,flag,d)
        else
          do j = 1,nrhs
            call ma64_solveD1(nfront,nelim,xlocal(:,j),flag,d)
          end do
        end if
      else
! data held in keep%rmain
        if (nrhs >= 4) then
          call ma64_solveD2(nfront,nelim,nrhs,xlocal,maxlen, &
               flag,keep%rmain(posdiag:posdiag+dreals-1))
        else
          do j = 1,nrhs
            call ma64_solveD1(nfront,nelim,xlocal(:,j),flag, &
                 keep%rmain(posdiag:posdiag+dreals-1))
          end do
        end if
      end if
      posdiag = posdiag + dreals

! copy xlocal back to x
      if (keep%index(1) >= 0) then
        do i = 1, nfront
          jvar = varlist(i)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
        end do
      else
        do i = 1, nfront
          jvar = keep%imain(loc+i-1)
          x(jvar,1:nrhs) = xlocal(i,1:nrhs)
        end do
      end if

      if (info%flag >= 0) return

480   info%iostat = keep%idata%iostat
      info%stat = keep%idata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return
    end subroutine diag_solve

!*******************************

! End solve

  end subroutine MA77_solve_double

!****************************************************************************
     subroutine MA77_enquire_posdef_double(d,keep,control,info)

! allows user to obtain the pivots used

 type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
 type (MA77_control), intent (inout) :: control ! See derived-type declaration
 type (MA77_info), intent (inout) :: info    ! See derived-type declaration

    real(wp), intent (out) :: d(:) ! d(i) will hold the the i-th pivot.

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!     for all path)
! -35 Call does not follow a call to MA77_factor with pos_def = .true.
! -37 Error in size of d

    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree
    integer(short)  :: depth ! depth in tree
    integer(short)  :: i
    integer(short)  :: ir ! recursive do loop variable
    integer(short)  :: flag ! local error flag
    integer(short)  :: kpiv
    integer(short)  :: n
    integer(short)  :: nchild ! number of children of node in tree
    integer(short)  :: nelt ! set to keep%nelt
    integer(short)  :: node ! node in tree
    integer(short)  :: nelim ! number of eliminations at node
    integer(short)  :: nout
    integer(short)  :: nroot ! number of roots (components)
    integer(short)  :: root ! root of tree
    integer(short)  :: pnode ! parent node in tree
    integer(long)   :: posdiag ! points to storage for diagonal entries
    integer(short)  :: st ! stat parameter

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_enquire_posdef'

    n     = keep%n
    nelt  = keep%nelt
! make checks
    if (keep%status /= 3) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (.not.keep%pos_def) then
      info%flag = -35
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (size(d) < n) then
      info%flag = -37
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    deallocate (count,stat=st)
    deallocate (cnode,stat=st)

    allocate (count(0:keep%maxdepth),cnode(0:keep%maxdepth),stat=st)
    if (st /= 0) go to 495

   nroot = size(keep%roots)
! initialise posdiag/locint
    posdiag = keep%posfac
    kpiv = 0
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. 
        call ask_posdef(node,kpiv)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 496
    end do ! roots


   go to 500

495   info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      go to 500

496   call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)

500 continue
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)

  contains

!*******************************

   subroutine ask_posdef(node,kpiv)

      integer(short), intent (in) :: node
      integer(short), intent (inout) :: kpiv

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

      nelim = keep%tree(node)%nelim
      if (nelim == 0) return

      if (keep%index(2) >= 0) then
! read required part of d
        call of01_read(keep%index(2),posdiag,nelim,d(kpiv+1:kpiv+nelim), &
             flag,keep%rdata,lp=-1)
        if (flag < 0) go to 485
!        write (6,*) d(kpiv+1:kpiv+nelim)
      else
! data held in keep%rmain
        d(kpiv+1:kpiv+nelim) = keep%rmain(posdiag:posdiag+nelim-1)
      end if
      do i = 1,nelim
        d(kpiv+i) = d(kpiv+i)**2
      end do
      posdiag = posdiag + nelim
      kpiv = kpiv + nelim

      if (info%flag >= 0) return

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

    end subroutine ask_posdef

    end subroutine ma77_enquire_posdef_double
!*************************************************
     subroutine MA77_enquire_indef_double(piv_order,d,keep,control,info)

! In indefinite case, the pivot sequence used will not necessarily be
! the same as that passed to ma77_factor (because of delayed pivots). This
! subroutine allows the user to obtain the pivot sequence that was
! actually used.
! also the entries of D^{-1} are returned using array d.

 type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
 type (MA77_control), intent (inout) :: control ! See derived-type declaration
 type (MA77_info), intent (inout) :: info    ! See derived-type declaration

    integer(short), intent (out) :: piv_order(:) !
!      If i is used to index a variable, its position in the pivot sequence will
!      be placed in piv_order(i), with its sign negative if it is
!      part of a 2 x 2 pivot; otherwise, piv_order(i) will be set to zero.
    real(wp), intent (out) :: d(:,:) ! The diagonal entries of D^{-1}
!      will be placed in d(1,:i) and the off-diagonal entries will 
!      be placed in d(2,:). The entries are held in pivot order.

! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!     for all path)
! -34 Call does not follow a call to MA77_factor with pos_def = .false.
! -36 Error in size of piv_order
! -37 Error in size of d

    integer(short) :: depth ! depth in tree
    integer(short) :: dreals ! set to 2*nelim
    integer(short) :: flag ! local error flag
    integer(short) :: i
    integer(short) :: ir ! recursive do loop variable
    integer(short) :: j ! temporary variable
    integer(short) :: jj ! temporary variable
    integer(short) :: kpiv
    integer(short) :: k ! temporary variable
    integer(short) :: maxlen ! holds largest number of variables in
!          an element list
    integer(short) :: n
    integer(short) :: nchild ! number of children of node in tree
    integer(short) :: nelt ! set to keep%nelt
    integer(short) :: node ! node in tree
    integer(short) :: nelim ! number of eliminations at node
    integer(short) :: nfront ! number of variables in node
    integer(short) :: nout
    integer(short) :: nroot ! number of roots (components)
    integer(short) :: pnode ! parent node in tree
    integer(short) :: root ! root of tree
    integer(short) :: st ! stat parameter
    integer(long)  :: loc ! location in file to read from
    integer(long)  :: locint
    integer(long)  :: posdiag ! points to storage for diagonal entries

    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree
    integer(short), allocatable :: varlist(:) ! used for variable lists
    real(wp), allocatable :: d_copy(:) ! used for reading in entries of D

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_enquire_indef'

    n     = keep%n
    nelt  = keep%nelt
! make checks
    if (keep%status /= 3) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (keep%pos_def) then
      info%flag = -34
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (size(piv_order) < n) then
      info%flag = -36
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if
    piv_order(1:n) = 0

    if (size(d,1) < 2 .or. size(d,2) < n) then
      info%flag = -37
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    deallocate (varlist,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (d_copy,stat=st)

    maxlen = keep%maxlen
    if (.not. keep%pos_def) maxlen = max(maxlen,keep%maxfrontb)
    allocate (count(0:keep%maxdepth),cnode(0:keep%maxdepth), &
      varlist(maxlen+1),d_copy(2*keep%maxfrontb),stat=st)
    if (st /= 0) go to 495

    nroot = size(keep%roots)
! initialise posdiag/locint
    posdiag = keep%posfac
    locint = keep%posint

    d(1:2,1:n) = zero ! ensure do not returned with this undefined

    kpiv = 0
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. 
        call ask_indef(node,kpiv)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 496
    end do ! roots

   go to 500

495   info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      go to 500

496   call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)

! Deallocate arrays we have done with.
500 deallocate (varlist,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (d_copy,stat=st)

  contains

!*******************************

   subroutine ask_indef(node,kpiv)

      integer(short), intent (in) :: node
      integer(short), intent (inout) :: kpiv

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

! The variable lists are held after the main integer lists
      nfront = keep%size_ind(node)
      loc = locint
! Increment locint so that it points to start of integer data for next node
      locint = locint + nfront + 1
      if (keep%index(1) >= 0) then
        call of01_read(keep%index(1),loc,nfront+1,varlist,flag, &
             keep%idata,lp=-1)
        if (flag < 0) go to 480
        nelim = varlist(nfront+1)
        if (nelim == 0) return
      else
! File not in use. Integers are held in keep%imain
        nelim = keep%imain(loc+nfront)
        if (nelim == 0) return
! Copy into varlist
        varlist(1:nelim) = keep%imain(loc:loc+nelim-1)
      end if

      dreals = 2*nelim
      if (keep%index(2) >= 0) then
! read required part of d_copy
        loc = posdiag
        call of01_read(keep%index(2),loc,dreals,d_copy, &
             flag,keep%rdata,lp=-1)
        if (flag < 0) go to 485
      else
! data held in keep%rmain
        d_copy(1:dreals) = keep%rmain(posdiag:posdiag+dreals-1)
      end if
      posdiag = posdiag + dreals

      j = 1
      jj = 1
      do i = 1,nelim
        if (j > nelim) exit
        kpiv = kpiv + 1
        k = varlist(j)
        if (d_copy(jj+1) == zero) then
! 1x1 pivot (which could be zero)
          d(1,kpiv) = d_copy(jj)
          piv_order(k) = kpiv
          j = j + 1; jj = jj + 2
        else
          d(1,kpiv)   = d_copy(jj)
          d(2,kpiv)   = d_copy(jj+1)
          d(1,kpiv+1) = d_copy(jj+2)
! negative flags for 2x2 pivots
          piv_order(k) = -kpiv
          kpiv = kpiv + 1
          k = varlist(j+1)
          piv_order(k) = -kpiv
          j = j + 2; jj = jj + 4
        end if
      end do

      if (info%flag >= 0) return

480   info%iostat = keep%idata%iostat
      info%stat = keep%idata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return
    end subroutine ask_indef

    end subroutine ma77_enquire_indef_double
!*************************************************
     subroutine MA77_alter_double(d,keep,control,info)

! In indefinite case, the entries of D^{-1} may be changed

 type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
 type (MA77_control), intent (inout) :: control ! See derived-type declaration
 type (MA77_info), intent (inout) :: info    ! See derived-type declaration

    real(wp), intent (in) :: d(:,:) ! The required diagonal entries
!      of D^{-1} must be placed in d(1,i) (i = 1,...n)
!      and the off-diagonal entries must be placed in d(2,i) (i = 1,...n-1).

! Possible error returns:
!  -1   Allocatation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!     for all path)
! -34 Call does not follow a call to MA77_factor with pos_def = .false.
! -37 Error in size of d

    integer(short)  :: depth ! depth in tree
    integer(short)  :: dreals ! set to 2*nelim
    integer(short)  :: flag ! local error flag
    integer(short)  :: i
    integer(short)  :: ir ! recursive do loop variable
    integer(short)  :: itemp(1)
    integer(short)  :: j ! temporary variable
    integer(short)  :: jj ! temporary variable
    integer(short)  :: kpiv
    integer(short)  :: n
    integer(short)  :: nchild ! number of children of node in tree
    integer(short)  :: nelt ! set to keep%nelt
    integer(short)  :: nfront ! number of variables in node
    integer(short)  :: node ! node in tree
    integer(short)  :: nelim ! number of eliminations at node
    integer(short)  :: nout
    integer(short)  :: nroot ! number of roots (components)
    integer(short)  :: pnode ! parent node in tree
    integer(short)  :: root ! root of tree
    integer(long)   :: loc ! location in file to read from
    integer(long)   :: locw ! location in file to write to
    integer(long)   :: locint
    integer(long)   :: posdiag ! points to storage for diagonal entries
    integer(short)  :: st ! stat parameter
    logical :: two ! set to .true. for 2x2 pivot

    integer(short), allocatable :: count(:) ! used for depth
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth
!      first search of tree
    real(wp), allocatable :: d_copy(:) ! used for writing
!      out changed entries on diagonal

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_enquire_indef'

    n     = keep%n
    nelt  = keep%nelt
! make checks
    if (keep%status /= 3) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (keep%pos_def) then
      info%flag = -34
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (size(d,1) < 2 .or. size(d,2) < n) then
      info%flag = -37
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (d_copy,stat=st)

    allocate (count(0:keep%maxdepth),cnode(0:keep%maxdepth), &
              d_copy(1:2*keep%maxfrontb),stat=st)
    if (st /= 0) go to 495

    nroot = size(keep%roots)
! initialise posdiag
    posdiag = keep%posfac
    locint = keep%posint

    kpiv = 0
    do ir = 1, nroot
      root = keep%roots(ir)
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
! Visit each node in depth-first search order
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. 
        call alter_indef(node,kpiv)
        if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        depth = depth - 1
        pnode = cnode(depth)
        node = pnode

      end do ! tree
      if (info%flag < 0) go to 496
    end do ! roots
   go to 500

495   info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      go to 500

496   call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)

! Deallocate arrays we have done with.
500 deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (d_copy,stat=st)

  contains

!*******************************
     subroutine alter_indef(node,kpiv)

      integer(short), intent (in) :: node
      integer(short), intent (inout) :: kpiv ! pivot count.

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)

! Immediate return for leaf node.
      if (node <= keep%nelt) return

! The number of eliminations actually performed at node is
! held in the main integer lists
      nfront = keep%size_ind(node)
      loc = locint + nfront
! Increment locint so that it points to start of integer data for next node
      locint = loc + 1
      if (keep%index(1) >= 0) then
        call of01_read(keep%index(1),loc,1,itemp,flag,keep%idata,lp=-1)
        if (flag < 0) go to 480
        nelim = itemp(1)
      else
! File not in use. Integers are held in keep%imain
        nelim = keep%imain(loc)
      end if
      if (nelim == 0) return

      j = 1
      jj = 1
      do i = 1,nelim
        if (j > nelim) exit
        kpiv = kpiv + 1
        two = .false.
! be careful! d(2,n) may not be defined
        if (kpiv < n) then
          if (d(2,kpiv) /= zero) two = .true.
        end if
        if (.not. two) then
! 1x1 pivot
          d_copy(jj)   = d(1,kpiv)
          d_copy(jj+1) = zero
          j = j + 1; jj = jj + 2
        else
          d_copy(jj)   = d(1,kpiv)
          d_copy(jj+1) = d(2,kpiv)
          d_copy(jj+2) = d(1,kpiv+1)
          d_copy(jj+3) = zero
          kpiv = kpiv + 1
          j = j + 2; jj = jj + 4
        end if
      end do

      dreals = 2*nelim
      if (keep%index(2) >= 0) then
! write required part of d
        locw = posdiag
        call of01_write(keep%index(2),locw,dreals,d_copy(1:dreals), &
             flag,keep%rdata,lp=-1)
        if (flag < 0) go to 485
      else
! data held in keep%rmain
        keep%rmain(posdiag:posdiag+dreals-1) = d_copy(1:dreals)
      end if
      posdiag = posdiag + dreals

      if (info%flag >= 0) return

480   info%iostat = keep%idata%iostat
      info%stat = keep%idata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

485   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16
      return

    end subroutine alter_indef

    end subroutine ma77_alter_double

!*************************************************

  subroutine MA77_scale_double(scale,keep,control,info,anorm)

! Scale matrix (equilibration).

! For details of keep, control, info : see derived type description
    type (MA77_keep), intent (inout) :: keep
    type (MA77_control), intent (in) :: control
    type (MA77_info), intent (inout) :: info
    real(wp), intent(out),optional :: anorm ! on exit, holds infinity norm of A

    real (wp), intent(out) :: scale(:) ! on exit, scale(1:n) holds the 
!     scaling factors.

! The following arrays are allocated by MA77_scale
    integer(short), allocatable :: ijw(:) ! 
    real (wp), allocatable :: dew(:) !
    real (wp), allocatable :: wnorm(:) ! Work array (only if anorm present)  
    real (wp), allocatable :: fa(:) ! frontal matrix
    real (wp), allocatable :: reals(:) ! used to hold a single col.
!            of the frontal matrix. Also used for reading in
!            a user-supplied element/row.
    integer(short), allocatable :: count(:) ! used for depth 
!      first search of tree to keep track of level we are at
    integer(short), allocatable :: cnode(:) ! used for depth 
!      first search of tree
    integer(short), allocatable :: varlist(:) ! used for holding
!            the variables at nodes in the tree.
!            Allocated to be of length keep%maxlen.
    integer(short), allocatable :: rlist(:) ! used for holding
!            the row and column indices at nodes in the tree.
!            They are allocated to be of length keep%maxlen. 
    integer(short), allocatable :: map(:) ! used to map between variables in
!            a child and those in its parent.
    integer(short), allocatable :: pos(:) ! Used for mapping into the front
!            (only array of length n)
!            Also used to hold column permutation in partial factorization
    integer(long), allocatable :: ip(:) ! Used to hold
!            positions of columns in the front.

  integer(short) :: cfront ! frontsize at node
  integer(short) :: cnvar ! number of (uneliminated) variables at child node
  integer(short) :: depth ! depth in tree
  integer(short) :: flag
  integer(short) :: ii ! temporary variable
  integer(short) :: i ! temporary variable
  integer(short) :: ic ! indicates which child of pnode node is
  integer(long)  :: ilong ! temporary variable
  integer(short) :: ir ! temporary variable
  integer(long)  :: ipp ! points to start of a column
  integer(long)  :: ip1 ! points to position after delayed cols in frontal mx
  integer(long)  :: ipc ! points to position of start of frontal mx for cnode
  integer(long)  :: ipj2 ! temporary variable (used assemble)
  integer(long)  :: ipk  ! temporary variable (used assemble)
  integer(short) :: ivar ! a variable
  integer(short) :: j1 ! temporary variable
  integer(short) :: j2 ! temporary variable
  integer(short) :: jc ! temporary variable
  integer(short) :: jcnode ! jc-th child node
  integer(short) :: jj ! temporary variable
  integer(short) :: j ! temporary variable
  integer(long)  :: jlong  ! temporary variable
  integer(short) :: jvar ! a variable
  integer(short) :: k  ! temporary variable
  integer(short) :: kk ! temporary variable
  integer(long)  :: klong ! temporary variable
  integer(short) :: lfa ! size of array fa
  integer(short) :: lreals ! length of array reals
  integer(short) :: l  ! temporary variable
  integer(long)  :: ln ! 
  integer(long)  :: length ! length of list to be read
  integer(long)  :: loc  ! location in real or integer superfile when reading
  integer(long)  :: locw ! location in real or integer superfile when writing
  integer(short) :: lw
  integer(short) :: maxlen ! set to keep%maxlen
  integer(short) :: maxit ! number of iterations
  integer(short) :: mvar ! copy of keep%mvar
  integer(short) :: nchild ! number of children of root
  integer(short) :: n ! set to keep%n
  integer(short) :: niter ! iterations of equilibration algorithm
  integer(short) :: nelim_jc
  integer(short) :: ncand ! number of variables that are candidates
!         for elimination at node
  integer(short) :: nfront ! order of frontal matrix
  integer(short) :: node ! node in tree
  integer(short) :: non_leaf ! number of children of node that are not
!         leaf nodes
  integer(short) :: nout ! output unit for errors
  integer(short) :: nroot ! number of roots (components) of the tree
  integer(short) :: nvar ! no. variables associated with node during analyse
  integer(long)  :: nvar_long ! no. variables associated with node 
!         during analyse
  integer(short) :: pnode ! parent node in tree
  integer(short) :: root ! root node
  integer(long)  :: size_fa ! size of frontal array fa
  integer(short) :: splitp ! split point at which space allocated for
!         generated frontal matrix
  integer(short) :: st ! stat parameter
  integer(long) :: loc1
  integer(long) :: rtop ! points to the top of the main real stack

  real(wp) :: err
  real(wp) :: thresh
  real (wp) :: multiplier 


! Possible error returns:
!  -1   Allocation error
!  -3   Error in sequence of calls
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -14   data for last element passed to MA77_input_reals was incomplete.
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open failed 
!       for all path)
! -22   For one or more elements, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.
! -30   front size is too large
! -39   size(scale) is too small


!   write(*,*) 'entering MA77_scale with keep%status=', keep%status

    keep%name = 'MA77_scale'

! Perform appropriate printing
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
      write (control%unit_diagnostics,'(//a)') &
       ' Entering MA77_scale'
      write (control%unit_diagnostics, &
     '(a,4(/a,i12),2(/a,i12),(/a,es12.4))') &
     ' control parameters (control%) :', &
     ' print_level         Level of diagnostic printing           = ', &
       control%print_level, &
     ' unit_diagnostics    Unit for diagnostics                   = ', &
       control%unit_diagnostics, &
     ' unit_error          Unit for errors                        = ', &
       control%unit_error, &
     ' unit_warning        Unit for warnings                      = ', &
       control%unit_warning, &
     ' maxit               Maximum number of iterations           = ', &
       control%maxit, &
     ' infnorm             Control for the scaling norm           = ', &
       control%infnorm, &
     ' thresh              Scaling stopping criteria              = ', &
       control%thresh

    end if

    nout = control%unit_error
    if (control%print_level < 0) nout = -1

    info%flag = 0

! Check status parameter
    if (keep%status < 2) then
      info%flag = -3
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    if (keep%n == 0) then
      if (present(anorm)) anorm = zero
      return
    end if

! check the last call to MA77_input_reals was complete.
    if (keep%inelrs > 0) then
! data for last element/row not complete
      info%flag = -14
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

! Deallocate arrays used by analyse but no longer needed.
    deallocate (keep%clist,stat=st)
    deallocate (keep%arow,stat=st)
    deallocate (keep%aelt,stat=st)
    deallocate (keep%iptr,stat=st)

    n = keep%n

    if (size(scale) < n) then
      info%flag = -39
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

    maxit = control%maxit
    if (maxit < 1) maxit = 1
    thresh = control%thresh

! Allocate arrays needed by scale
    multiplier = one
    mvar = keep%mvar
    maxlen = keep%maxlen

    deallocate (fa,stat=st)
    lfa = keep%maxfront
    lfa = (lfa*(lfa+1))/2
    if (lfa > keep%lup) then
      info%flag = -30
      call MA77_print_iflag(keep,nout,info%flag)
      go to 25
    end if

    allocate (fa(lfa),stat=st)
    if (st /= 0) then
! if arrays are being used, we try and switch to using files
      if (all(keep%index(1:4) > 0)) then
! files are being used so we have failed
        go to 20
      else
! at least one array is being used.
! Try and write the integers to file
        if (keep%index(1) < 0) then
          keep%index(1) = -keep%index(1)
          loc1 = 1
          lw = keep%ifree - 1
          call of01_write(keep%index(1),loc1,lw,keep%imain, &
             flag,keep%idata,lp=-1)
          if (flag < 0) then
            info%iostat = keep%idata%iostat
            info%stat = keep%idata%stat
            info%flag = flag
            if (flag == -17) info%flag = -16
            go to 25
          end if
          keep%used = keep%used - size(keep%imain)
          deallocate (keep%imain,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Try switching reals to a file
        if (keep%index(2) < 0) then
          keep%index(2) = -keep%index(2)
          loc1 = 1
          lw = keep%posfac - 1
          call of01_write(keep%index(2),loc1,lw,keep%rmain, &
             flag,keep%rdata,lp=-1)
          if (flag < 0) then
            info%iostat = keep%rdata%iostat
            info%stat = keep%rdata%stat
            info%flag = flag
            if (flag == -17) info%flag = -16
            go to 25
          end if
          keep%used = keep%used - size(keep%rmain)*2
          deallocate (keep%rmain,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Nothing yet written to the in-core arrays index(4) so switch to files
        if (keep%index(4) < 0) then
          deallocate(keep%rwdelay,stat=st)
          keep%index(4) = -keep%index(4)
          keep%used = keep%used - 2*control%storage_indef
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Now try to write stack to file (at moment, nothing in the stack)
        if (keep%index(3) < 0) then
          keep%index(3) = -keep%index(3)
          keep%used = keep%used - size(keep%rwork)*2
          deallocate (keep%rwork,stat=st)
          allocate (fa(lfa),stat=st)
          if (st == 0) go to 10
        end if
! Have not been able to allocate front
        info%flag = -30; info%stat = st
        call MA77_print_iflag(keep,nout,info%flag)
        go to 25
      end if
    end if
! frontal matrix has been allocated
 10 continue
    size_fa = lfa

    deallocate (ip,stat=st)
    allocate (ip(maxlen),stat=st)
    if (st /= 0) go to 20

    lreals = max(mvar,maxlen)

    deallocate (rlist,stat=st)
    deallocate (varlist,stat=st)
    deallocate (reals,stat=st)
    deallocate (pos,stat=st)
    deallocate (map,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)

    allocate (count(0:keep%maxdepth),cnode(0:keep%maxdepth),stat=st)
    if (st /= 0) goto 20

    allocate (reals(lreals),rlist(maxlen+1),varlist(maxlen), &
              map(maxlen),pos(1:n),stat=st)
    if (st /= 0) go to 20
! Initialise the array pos to zero
    pos = 0

! reset keep%rfree and keep%ifree
    keep%rfree = keep%posfac
    keep%ifree = keep%posint

    keep%rtopmx = 1_long

! Perform scaling.

    deallocate (ijw,stat=st)
    deallocate (dew,stat=st)
    allocate (ijw(n),dew(n),stat=st)
    if (st /= 0) go to 20

    niter = 0
    ijw(1:n) = 0
    scale(1:n) = one
    info%niter = niter

    if (present(anorm)) then
      deallocate (wnorm,stat=st)
      allocate (wnorm(n),stat=st)
      if (st /= 0) go to 20
      wnorm(1:n) = zero
    end if

iter: do

    niter = niter + 1

    dew(1:n) = zero

    nroot = size(keep%roots)
    do ir = 1, nroot
      root = keep%roots(ir)
! Initialise the top of the real stack
      rtop = 1
! Visit each node in depth-first search order
      node = root
      depth = 1
      count(depth) = 0 ! no. of visited children
      do while(depth.gt.0)
         cnode(depth) = node
! Visit children
         nchild = size(keep%tree(node)%child)
         if (count(depth) < nchild) then
            count(depth) = count(depth) + 1
            if (keep%tree(node)%child(count(depth)) > keep%nelt) then
! The next child of node is not a leaf node so descend to it
               node = keep%tree(node)%child(count(depth))
               depth = depth + 1
               count(depth) = 0
               cycle
            end if
         end if

! Perform work at node. Assemble contributions
! from children of node that are leaf nodes 
! and then perform the work at node.
! Then either stack the contribution from node or
! assemble it into its parent (node is the ic-th child of its parent)
         if (niter == 1 .and. present(anorm)) then
           call scale77(node,n,cfront,scale,dew,ijw,wnorm)
         else
           call scale77(node,n,cfront,scale,dew,ijw)
         end if
         if (info%flag < 0 .or. depth <= 1) exit
         depth = depth - 1
         pnode = cnode(depth) 
         ic = count(depth)
         call assemble(pnode,node,ic,rtop)
         if (info%flag < 0) exit
! If we have reached the root, we are done
        if (node == root) exit

! Go to parent
        node = pnode

      end do ! tree

      if (info%flag < 0) then
        keep%status = -1
        call MA77_print_iflag(keep,nout,info%flag,st=info%stat, &
             ios=info%iostat)
        go to 25
      end if
    end do

       if (niter == 1 .and. present(anorm)) then
         anorm = wnorm(1)
         do i = 2,n
           anorm = max(anorm,wnorm(i))
         end do
       end if

! Accumulate scaling factors
       if (control%infnorm == 0) then
! infinity norm
         do i = 1,n
           if (ijw(i) > 0) scale(i) = scale(i)*sqrt(dew(i))
         end do

         k = 0
! Loop over columns, checking to see if maximum in column j
! and row j is on the diagonal
         do j = 1,n
! ijw(j) holds row index of largest entry in col. j
           i = ijw(j)
           if (i > 0) then
! ijw(i) holds column index of largest entry on row i
             if (ijw(i) == j) then
               ijw(i) = -ijw(i)
               k = k + 1
               if (i /= j) then
                 ijw(j) = -ijw(j)
                 k = k + 1
               end if
             end if
           else
             k = k + 1
           end if
         end do
         if (k == n) exit
       else
! 1- norm
         do i = 1,n
           if (ijw(i) == 0) cycle
           scale(i) = scale(i)*sqrt(dew(i))
         end do

      end if

      err = zero
      do i = 1,n
        if (ijw(i) > 0) err = max( err, abs(one-dew(i)) )
      end do

      if ( err <= thresh ) exit
      if (niter == maxit) exit

    end do iter
    info%niter = niter

! Store the inverse scaling factor
    do i = 1,n
      if (ijw(i) == 0) cycle
      scale(i) = one/scale(i)
    end do 

20  if (st /= 0) then
      info%flag = -1; keep%status = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
    end if

 25 deallocate (reals,stat=st)
    deallocate (map,stat=st)
    deallocate (pos,stat=st)
    deallocate (rlist,stat=st)
    deallocate (varlist,stat=st)
    deallocate (fa,stat=st)
    deallocate (count,stat=st)
    deallocate (cnode,stat=st)
    deallocate (ijw,stat=st)
    deallocate (dew,stat=st)
    deallocate (wnorm,stat=st)
    deallocate (ip,stat=st)

    info%index(1:4) = keep%index(1:4)

    keep%flag = info%flag
    if (info%flag < 0) return

    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) then
        write (control%unit_diagnostics,'(/a)') &
       ' Completed scaling'
    end if
    if (present(anorm)) then
      if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
        write (control%unit_diagnostics,'(/a,es12.4)') &
      ' Leaving MA77_scale with        anorm                    = ', anorm
    end if

contains

!*******************************

   subroutine scale77(node,n,cfront,de,dew,ijw,wnorm)

! This subroutine assembles contributions from any children of node
! that are leaf nodes, and then finds largest entries in fully summed
! part, and optionally computes row sums for fully summed rows.

    integer(short), intent(in) :: node ! node in tree
    integer(short), intent(in) :: n
    integer(short), intent(out) :: cfront ! frontsize at node
    integer(short), intent(inout) :: ijw(n)
    real(wp), intent(in) :: de(n) ! accumulated scaling factors
    real(wp), intent(inout) :: dew(n) ! scaling factors on current iteration
    real(wp), optional, intent(inout) :: wnorm(n) ! work array for
!             computing matrix norm

! Possible error returns:
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -15   Error from of01 (error in Fortran write)
! -22   For one or more elements/rows, MA77_input_vars was called but
!       no corresponding call was made to MA77_input_reals.

      nvar = abs(keep%size(node))
      if (nvar == 0) return
! If node is a leaf, immediate return.
      if (node <= keep%nelt) return

      nvar_long = nvar
      nfront = nvar
      ln = (nvar_long*(nvar_long+1))/2
      ip1 = 1 + size_fa - ln
      lfa = ln

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild
! In element case, add up the number of children of node that are leaf nodes
! .... they are the last children (in row case, leaf nodes are original
! rows and they are not part of the tree)
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if
! If all the children of node are leaf children, initialize front
      if (non_leaf == 0) fa(ip1:ip1+lfa-1) = zero

! Read the variables for node from the main integer superfile into 
! rlist and set mapping into front 
      loc = keep%ifile(node)
      call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
           rlist(1:nvar),flag,keep%idata,-1)
      if (flag < 0) go to 480
      do i = 1,nvar
        rlist(i) = abs(rlist(i))
        ivar = rlist(i)
        pos(ivar) = i
      end do

! In element case, must now deal with children that are leaf nodes
! and in row case, deal with original rows

! Set ip(j) to point to start of jth col.
      ip(1) = ip1
      do j = 2,nvar
        ip(j) = ip(j-1) + nvar - j + 2
      end do

      if (keep%element_input) then
        do ic = non_leaf+1,nchild
          jcnode = keep%tree(node)%child(ic)
! Check reals were entered.
          if (keep%rfile(jcnode) == -1 .and. keep%size(jcnode) /= 0) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
          cnvar = abs(keep%size(jcnode))
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
               varlist,flag,keep%idata,-1)
          if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos.
          do j = 1,cnvar
            k = varlist(j)
            map(j) = pos(k)
          end do

! Read reals from main file one column at a time
          loc = cnvar
          loc = keep%rfile(jcnode) - (loc*loc+loc)/2
          do j = 1,cnvar
            lreals = cnvar - j + 1
            call MA77_read_real(keep%index(2),keep%rmain,loc,lreals, &
                 reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480

! Add the reals from the original element into the front
            k = 1
            jj = map(j)
            do i = j, cnvar
              ii = map(i)
              if (ii >= jj) then
                klong = ip(jj) + ii - jj
              else
                klong = ip(ii) + jj - ii
              end if
              fa(klong) = fa(klong) + reals(k)
              k = k + 1
            end do
            loc = loc + lreals
          end do
        end do

      else
! For row entry, add into the front the original rows corresponding
! to the eliminated variables.
        do j = 1, keep%tree(node)%nelim
          jcnode = rlist(j)
          cnvar = abs(keep%size(jcnode))
          if (cnvar == 0) cycle
          if (keep%rfile(jcnode) == -1) then
! User has not entered required reals
            info%flag = -22; return
          end if
          loc = keep%ifile(jcnode)
! Read in the variable list for row jcnode and its reals
          call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar,varlist,&
               flag,keep%idata,-1)
          if (flag < 0) go to 480
! Recall that keep%rfile(jcnode) points to position in the main reals
! file after the last entry in the list of reals for jcnode
          loc = keep%rfile(jcnode) - cnvar
          call MA77_read_real(keep%index(2),keep%rmain,loc,cnvar,reals,flag,&
               keep%rdata,-1)
          if (flag < 0) go to 480

! Loop over entries in the row jcnode, putting
! them into row/col j of the frontal matrix (lower triangular part).
! ipp points to start of col. j
           jlong = j
           ipp = ip1 + (nvar_long+1)*(jlong-1) - (jlong*(jlong-1))/2
           do jj = 1, cnvar
             jvar = varlist(jj)
             i = pos(jvar)
             if (i < j) cycle
! Put entry into col j, row i.
             fa(ipp+i-j) = fa(ipp+i-j) + reals(jj)
           end do
        end do
      end if

! Restore pos to zero
      do i = 1,nvar
        ivar = rlist(i)
        pos(ivar) = 0
      end do

! Have now assembled the frontal matrix and are ready to find
! largest entries in fully summed rows/cols. It is of order nfront
! but we hold only the lower triangular part.
! ncand is no. of fully summed rows/cols

     ncand = keep%tree(node)%nelim

     if (control%infnorm == 0) then
! infinity norm
       call getmx(nfront,ncand,fa(ip1),rlist,n,de,dew,ijw)
     else
       call get_norm1(nfront,ncand,fa(ip1:ip1),rlist,n,de,dew,ijw)
     end if
     if (present(wnorm)) &
       call get_anorm(nfront,ncand,fa(ip1),rlist,n,wnorm)

     cfront = nfront

     if (info%flag >= 0) return

480  info%iostat = keep%rdata%iostat
     info%stat = keep%rdata%stat
     info%flag = flag

     return

    end subroutine scale77

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine assemble(node,cnode,ic,rtop)

! This subroutine either stacks the contribution from cnode
! or assembles it into its parent node. cnode is the ic-th child of node.
! (it is same as assemble54)

    integer(short), intent(in) :: node ! node in tree
    integer(short), intent(in) :: cnode ! child of node
    integer(short), intent(in) :: ic ! position of ic in sibling list
    integer(long), intent(inout) :: rtop ! points to top of main real stack

! Possible error returns:
!  -1   Allocation error
!  -5   Error from of01 (error in Fortran inquire)
!  -6   Error from of01 (error in Fortran read)
!  -7   Error from of01 (error in Fortran open)
! -15   Error from of01 (error in Fortran write)
! -16   -17 Error from of01_read or of01_write (Fortran open 
!       failed for all path)

      nchild = size(keep%tree(node)%child)
      non_leaf = nchild
! In element case, add up the number of children of node that are leaf nodes
      if (keep%element_input) then
        if (nchild == 0) return
        do i = nchild,1,-1
          jcnode = keep%tree(node)%child(i)
          if (jcnode > keep%nelt) exit
        end do
        non_leaf = i
      end if

! set ipc to point to first location in fa that was used for cnode
        cnvar = keep%size(cnode) - keep%tree(cnode)%nelim
        length = cnvar
        length = (length*length + length)/2
        ipc = 1 + size_fa - length

        nvar = abs(keep%size(node))
        nfront = nvar
        nvar_long = nvar
        ln = (nvar_long*nvar_long + nvar_long)/2
! set ip1 to point to first location in fa that will be used for front at node
        ip1 = 1 + size_fa - ln
        lfa = ln
        splitp = keep%splitp(node)
! If we are not yet at the split point, write generated element for cnode onto
! top of the real stack and return
        if (ic < splitp) then
          locw = rtop
          lreals = cnvar
          ipp = ipc
          do j = 1,cnvar
            call MA77_write_real(keep%index(3),keep%size_rwork,keep%rwork, &
                 locw,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,        &
                 keep%maxstore,keep%used)
            if (flag < 0) go to 480
            ipp = ipp + lreals
            locw = locw + lreals
            lreals = lreals - 1
          end do
! Adjust top of stack
          rtop = rtop + length
          keep%rtopmx = max(keep%rtopmx,rtop)
          return
        end if

! Read the variables for node from the main integer superfile.
        loc = keep%ifile(node)
        call MA77_read_integer(keep%index(1),keep%imain,loc,nvar,&
             rlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
        do i = 1,nvar
          ivar = abs(rlist(i))
          pos(ivar) = i
        end do

! Set ip(j) to point to start of col j in array fa.
        ip(1) = ip1
        do j = 2,nvar
          ip(j) = ip(j-1) + nvar - j + 2
        end do

! read integer list for cnode from main integer file
        loc = keep%ifile(cnode) + keep%tree(cnode)%nelim
        call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
             varlist,flag,keep%idata,-1)
        if (flag < 0) go to 480
! Set mapping from cnode into the front
        do j = 1,cnvar
          k = abs(varlist(j))
          map(j) = pos(k)
        end do

        if (ic == splitp .or. ic == non_leaf) then
          if (cnvar /= nvar) then
! Expand cnode into the front
            j1 = 1
            do j = 1,cnvar
              j2 = map(j)
! col. j of cnode maps to col. j2 of front
! set cols j1:j2-1 of front to zero
              fa(ip(j1):ip(j2)-1) = zero
! Expand col. j of cnode into col. j2 of front.
! col j of cnode is currently at col k=j+nvar-cnvar of front
              k = j + nvar - cnvar
              if (j2 /= k) then
                ipj2 = ip(j2)
                ipk = ip(k)
                fa(ipj2:ip(j2+1)-1) = zero
                do l = j,cnvar
                  i = map(l)
                  fa(ipj2+i-j2) = fa(ipk+l-j)
                end do
              else
! j2 = k, so all columns are now in place
                j1 = nvar + 1
                exit
              end if
              j1 = j2 + 1
            end do
! Check final cols are set to 0.
            klong = size_fa
            if (j1 <= nvar) fa(ip(j1):klong) = zero
          end if
        end if

        if (ic == splitp) then
! At split point.
! Merge children 1:splitp-1 into rows/columns 1:nfront of frontal matrix
          do 160 jc = splitp-1, 1, -1
            jcnode = keep%tree(node)%child(jc)
            nelim_jc = keep%tree(jcnode)%nelim
            loc = keep%ifile(jcnode) + nelim_jc
            cnvar = keep%size(jcnode) - nelim_jc
            call MA77_read_integer(keep%index(1),keep%imain,loc,cnvar, &
                 varlist,flag,keep%idata,-1)
            if (flag < 0) go to 480
! Set mapping from jcnode into the front using array pos (length n)
            do j = 1,cnvar
              k = abs(varlist(j))
              map(j) = pos(k)
            end do
! Read in reals for child from the top of stack, one col. at a time,
! using map to read directly into the front.
            length = cnvar
            length = (length*length + length)/2
            loc = rtop - length
            lreals = cnvar
! move pointer to top of stack
            rtop = loc
            do j = 1, cnvar
              ilong = map(j)
              ipp = ip1 + nvar_long*(ilong-1) - (ilong*(ilong-1))/2
              call MA77_read_discard_real(keep%index(3),keep%rwork,loc,lreals,&
                   fa(ipp:ipp+nvar-1),flag,keep%rdata,-1,map=map(j:cnvar))
              if (flag < 0) go to 480
              loc = loc + lreals
              lreals = lreals - 1
            end do
! End of loop over children
160       continue

! At this point fa holds the reals in the frontal matrix in packed form.
! If split point is not also last non-leaf child, write fa onto the stack
          if (ic /= non_leaf) then
            lreals = nvar
            loc = rtop
            ipp = ip1
            do j = 1,nvar
              call MA77_write_real(keep%index(3),keep%size_rmain,keep%rwork,   &
                   loc,lreals,fa(ipp:ipp+lreals-1),flag,keep%rdata,-1,         &
                   keep%maxstore,keep%used)
              if (flag < 0) go to 480
              ipp = ipp + lreals
              loc = loc + lreals
              lreals = lreals - 1
            end do
            rtop = rtop + (nvar_long*nvar_long+nvar_long)/2
            keep%rtopmx = max(keep%rtopmx,rtop)
          end if

        else if (ic < non_leaf) then
! We have passed split point and not yet reached last non-leaf child.
! The reals for jcnode are still in fa.
! Loop over the cols of cnode, merging the child into its parent.
            do j = 1,cnvar
              i = map(j)
! Read column i of the frontal matrix from the stack. Find its position.
               loc = rtop - (nvar*nvar+nvar)/2
! loc now points to the first entry in the front on the stack.
! Set loc to point to first entry in col. i and
! lreals to number of entries in col. i.
              loc = loc + nvar*(i-1) - ((i-2)*(i-1))/2
              lreals = nvar - i + 1
              call MA77_read_real(keep%index(3),keep%rwork,loc,lreals,reals,&
                   flag,keep%rdata,-1)
              if (flag < 0) go to 480
! ipp points to start of col. j of jcnode, which is held in fa
              jlong = j
              ipp = ipc + cnvar*(jlong-1) - ((jlong-2)*(jlong-1))/2
              do k = j,cnvar
                kk = map(k)
                reals(1+kk-i) = reals(1+kk-i) + fa(ipp+k-j)
              end do
! Write the updated reals back onto the stack
              call MA77_write_real(keep%index(3),keep%size_rmain,keep%rwork, &
                  loc,lreals,reals,flag,keep%rdata,-1,keep%maxstore,keep%used)
              if (flag < 0) go to 480
            end do

        else if (ic == non_leaf) then
! Dealing with last non-leaf child (and it is not split point)
! Read frontal matrix from top of stack, one column at a time and add
! in with last non-leaf child (which we have already expanded in fa).
          loc = rtop - (nvar_long*nvar_long+nvar_long)/2
          rtop = loc
          ipp = ip1
          lreals = nvar
          do j = 1,nvar
            call MA77_read_discard_real(keep%index(3),keep%rwork,loc,&
                 lreals,reals,flag,keep%rdata,-1)
            if (flag < 0) go to 480
            fa(ipp:ipp+lreals-1) = fa(ipp:ipp+lreals-1) + reals(1:lreals)
            ipp = ipp + lreals
            loc = loc + lreals
            lreals = lreals - 1
          end do

        end if

! reset pos
      do i = 1,nvar
        ivar = abs(rlist(i))
        pos(ivar) = 0
      end do

      if (info%flag >= 0) return

480   info%iostat = keep%rdata%iostat
      info%stat = keep%rdata%stat
      info%flag = flag
      if (flag == -17) info%flag = -16 

    end subroutine assemble
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getmx(m,p,a,rlist,n,de,dew,ijw)

! Accumulate the maximum entry in each row and column.
! Entries in a are unscaled
! de holds accumulated scaling factors (initialised to 1 before first call)

      integer(short), intent(in) :: m ! order frontal matrix
      integer(short), intent(in) :: n ! global order
      integer(short), intent(in) :: p ! number fully summed rows/cols
      integer(short), intent(in) :: rlist(m) ! global indices
      integer(short), intent(inout) :: ijw(n)
      real(wp), intent(in) :: a(*)
      real(wp), intent(in) :: de(n) ! accumulated scaling factors
      real(wp), intent(inout) :: dew(n) ! scaling factors on current iteration

      integer :: i,j,k
      integer :: ivar,jvar
      real(wp) :: dtemp,s

! Loop over first p columns
        k = 1
        do j = 1,p
          jvar = rlist(j)
          dtemp = de(jvar)
          if (ijw(jvar) >= 0) then
! Find maximum entry in the column on or below diag.
            do i = j,m
              ivar = rlist(i)
              s = abs(a(k)) / (de(ivar)*dtemp)
              if (dew(jvar) < s) then
                dew(jvar) = s
                ijw(jvar) = ivar
              end if
              if (ijw(ivar) >= 0) then
                if (dew(ivar) < s) then
                  dew(ivar) = s
                  ijw(ivar) = jvar
                end if
              end if
              k = k + 1
            end do
          else
            do i = j,m
              ivar = rlist(i)
              if (ijw(ivar) >= 0) then
                s = abs(a(k)) / (de(ivar)*dtemp)
                if (dew(ivar) < s) then
                  dew(ivar) = s
                  ijw(ivar) = jvar
                end if
              end if
              k = k + 1
            end do
          end if
        end do

      end subroutine getmx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_norm1(m,p,a,rlist,n,de,dew,ijw)

! Accummulate 1-norm for the first p rows/columns

! Entries in a are unscaled
! de holds accumulated scaling factors (initialised to 1 before first call)

      integer(short), intent(in) :: m ! order frontal matrix 
      integer(short), intent(in) :: n ! global order
      integer(short), intent(in) :: p ! number fully summed rows/cols
      integer(short), intent(in) :: rlist(m) ! global indices
      integer(short), intent(inout) :: ijw(n)
      real(wp), intent(in) :: a(*)
      real(wp), intent(in) :: de(n) ! accumulated scaling factors
      real(wp), intent(inout) :: dew(n) ! scaling factors on current iteration

      integer :: i,j,k
      integer :: ivar,jvar
      real(wp) :: s

      do i = 1,m
        ivar = rlist(i)
! flag ivar as having been used
        ijw(ivar) = 1
      end do

        k = 1
        do j = 1,p
          jvar = rlist(j)
! deal with diagonal entry
          s = abs(a(k)) / (de(jvar)*de(jvar))
          dew(jvar) = dew(jvar) + s
          k = k + 1
! deal with entries below diagonal
          do i = j+1,m
            ivar = rlist(i)
            s = abs(a(k)) / (de(ivar)*de(jvar))
            dew(jvar) = dew(jvar) + s
            dew(ivar) = dew(ivar) + s
            k = k + 1
          end do
        end do

      end subroutine get_norm1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_anorm(m,p,a,rlist,n,wnorm)

! Accumulate row sums

      integer(short), intent(in) :: m ! order frontal matrix 
      integer(short), intent(in) :: n ! global order
      integer(short), intent(in) :: p ! number fully summed rows/cols
      integer(short), intent(in) :: rlist(m) ! global indices
      real(wp), intent(in) :: a(*)
      real(wp), intent(inout) :: wnorm(n)

      integer :: i,j,k
      integer :: ivar,jvar

! Loop over first p cols
        k = 1
        do j = 1,p
          jvar = rlist(j)
          wnorm(jvar) = wnorm(jvar) + abs(a(k))
          k = k + 1 
          do i = j+1,m
            wnorm(jvar) = wnorm(jvar) + abs(a(k))
            ivar = rlist(i)
            wnorm(ivar) = wnorm(ivar) + abs(a(k))
            k = k + 1   
          end do
        end do

      end subroutine get_anorm



! End scale

  end subroutine MA77_scale_double

!*************************************************

  subroutine MA77_print_iflag(keep,nout,iflag,ie,st,ios)
    integer(short), intent (in) :: iflag, nout
    integer(short), intent (in), optional :: ie, st, ios

    type (MA77_keep), intent (in) :: keep

    if (nout < 0) return
    if (iflag < 0) then
      write (nout,'(/3a,i3)') ' Error return from ',trim(keep%name),&
                             '. Error flag = ', iflag
      if (present(ie)) write (nout,'(a,i8)') &
        ' Error on input of element/row = ', ie
    else
      write (nout,'(/3a,i3)') ' Warning from ',trim(keep%name),&
         '. Warning flag = ', iflag
      if (present(ie)) write (nout,'(a,i8)') &
        ' Warning on input of element/row = ', ie
    end if

! Errors
    if (iflag == -1) then
      if (present(st)) write (nout,'(a,i3)') &
        ' Allocation error. stat parameter = ', st
    else if (iflag == -3) then
      write (nout,'(a)') ' Error in sequence of calls.'
    else if (iflag == -4) then
      write (nout,'(a,i10)') ' n = ', keep%n
    else if (iflag == -5) then
      write (nout,'(a)') ' Error in Fortran inquire'
      if (present(st)) write (nout,'(a,i3)') &
        ' iostat parameter = ', ios
    else if (iflag == -6) then
      write (nout,'(a)') ' Error in Fortran read'
      if (present(st)) write (nout,'(a,i3)') &
        ' iostat parameter = ', ios
    else if (iflag == -7) then
      write (nout,'(a)') ' Error in Fortran open'
      if (present(st)) write (nout,'(a,i3)') &
        ' iostat parameter = ', ios
    else if (iflag == -9) then
      write (nout,'(a,i10)') &
    ' MA77_input_reals has already been called for element/row = ', ie
    else if (iflag == -11) then
      write (nout,'(a)') ' Either pos-def = .true. but matrix not positive or'
      write (nout,'(a)') ' pos-def = .false. and matrix found to be singular'
    else if (iflag == -12) then
!      write (nout,'(a)') ' A file with a name supplied by user already exists'
    else if (iflag == -13) then
      write (nout,'(a)') ' len(filename) (or len(restart_file)) is too large'
    else if (iflag == -14) then
      write (nout,'(a)') &
        ' Data for last element/row input to MA77_input_reals was incomplete'
    else if (iflag == -15) then
      write (nout,'(a)') ' Error in Fortran write'
      if (present(st)) write (nout,'(a,i3)') &
        ' iostat parameter = ', ios
    else if (iflag == -16) then
      write (nout,'(a)') &
 ' len(path) too large or Fortran open unsuccessful for all elements of path'
    else if (iflag == -17) then
      write (nout,'(a)') &
    ' Duplicates and/or out-of-range entries were detected by MA77_input_vars.'
      write (nout,'(a)') &
    ' Thus all reals for each element/row must be entered in a single call.'
    else if (iflag == -18) then
      write (nout,'(a,i10)') ' nelt = ', keep%nelt
    else if (iflag == -19) then
      write (nout,'(a)') ' length <  0'
    else if (iflag == -20) then
      write (nout,'(a,i10)') ' job out of range.'
    else if (iflag == -21) then
      write (nout,'(a)') ' Error in user-supplied elimination order'
    else if (iflag == -22) then
      write (nout,'(a)') &
        ' MA77_input_reals was not called for one or more element/row'
    else if (iflag == -23) then
      write (nout,'(a)') ' control%buffer_lpage(1:2) out of range'
    else if (iflag == -24) then
      write (nout,'(a)') ' Error in size of x'
    else if (iflag == -25) then
      write (nout,'(a)') ' Error in size of resid'
    else if (iflag == -26) then
      write (nout,'(a)') ' control%buffer_npage(1:2) < 1'
    else if (iflag == -27) then
      write (nout,'(a)') ' control%maxstore out of range'
    else if (iflag == -28) then
      write (nout,'(a)') ' A file that is required does not exist'
    else if (iflag == -29) then
      write (nout,'(a)') &
    ' pos_def = .true. but 2x2 pivots were supplied on the call to MA77_analyse'
    else if (iflag == -30) then
      write (nout,'(a)') ' The front size is too large for allocation'
    else if (iflag == -31) then
      write (nout,'(a)') ' All entries in element/row out-of-range'
    else if (iflag == -32) then
      write (nout,'(a)') ' Too many reals input for this element/row'
    else if (iflag == -33) then
      write (nout,'(a)') ' Number of entries in variable list is < 0'
    else if (iflag == -34) then
      write (nout,'(a)') ' On call to MA77_factor, pos_def = .true.'
    else if (iflag == -35) then
      write (nout,'(a)') ' On call to MA77_factor, pos_def = .false.'
    else if (iflag == -36) then
      write (nout,'(a)') ' Error in size of piv_order'
    else if (iflag == -37) then
      write (nout,'(a)') ' Error in size of d'
    else if (iflag == -38) then
      write (nout,'(a)') &
    ' control%static < control%small and control%static /= 0.0'
    else if (iflag == -39) then
      write (nout,'(a)') ' Error involving scale'
    else if (iflag == -40) then
      write (nout,'(a)') ' IEEE infinities detected. Factorization terminated.'
    else if (iflag == -41) then
      write (nout,'(a)') ' Error in Fortran CLOSE statement.'


! Warnings
    else if (iflag == 4) then
      write (nout,'(a)') ' Matrix found to be singular'
    end if

  end subroutine MA77_print_iflag

!************************************************************************
  subroutine MA77_read_integer(ifile,keep_array,loc,length,read_array, &
             flag,data,lp)

! If files are in use (ifile >= 0), call of01_read to read into read_array.
! otherwise, copy from keep_array into read_array.

   integer(short), intent(in) :: ifile ! unit for reading from
   integer(short), allocatable :: keep_array(:) ! if files not in use,
!          we copy from this array into read_array
   integer(long), intent(in) :: loc ! position in keep_array to read from
   integer(short), intent(in) :: length ! number of reals to be read
   integer(short), dimension(:), intent(out) :: read_array
!          read into this array
   integer(short), intent(out) :: flag ! error flag
   type (of01_idata), intent(inout) :: data
   integer(short), intent(in) :: lp ! unit for of01 error messages

   integer(short) :: i,k

   flag = 0
   if (length == 0) return

   if (ifile >= 0) then
! read using of01
     call of01_read(ifile,loc,length,read_array,flag,data,lp=lp)
   else
! Copy from keep_array into read_array
! use loop unrolling to speed up the copying
     k = mod(length,7)
     read_array(1:k) = keep_array(loc:loc+k-1)
     do i = k+1,length, 7
       read_array(i) = keep_array(loc-1+i)
       read_array(i+1) = keep_array(loc+i)
       read_array(i+2) = keep_array(loc+i+1)
       read_array(i+3) = keep_array(loc+i+2)
       read_array(i+4) = keep_array(loc+i+3)
       read_array(i+5) = keep_array(loc+i+4)
       read_array(i+6) = keep_array(loc+i+5)
     end do
!    read_array(1:length) = keep_array(loc:loc+length-1)
   end if

  end subroutine MA77_read_integer

!************************************************************************

  subroutine MA77_read_real(rfile,keep_array,loc,length,read_array, &
             flag,data,lp,map)

! If files are in use (rfile >= 0), call of01_read to read into read_array.
! otherwise, copy from keep_array into read_array.

   integer(short), intent(in) :: rfile ! unit for reading from
   real(wp), allocatable :: keep_array(:) ! if files not in use,
!       we copy from this array into read_array
   integer(long), intent(in) :: loc ! position in keep_array to read from
   integer(short), intent(in) :: length ! number of reals to be read
   real(wp), dimension(:), intent(inout) :: read_array ! read into this array
   integer(short), intent(out) :: flag ! error flag
   type (of01_rdata), intent(inout) :: data
   integer(short), intent(in) :: lp ! unit for of01 error messages
   integer(short), optional, intent (in) :: map(length)

! Local scalars
   integer(short) :: i,j,k

   flag = 0
   if (length == 0) return

   if (rfile >= 0) then
! read using of01
     if (present(map)) then
       call of01_read(rfile,loc,length,read_array,flag,data,lp=lp,map=map)
     else
       call of01_read(rfile,loc,length,read_array,flag,data,lp=lp)
     end if
   else
! Copy from keep_array into read_array
     if (present(map)) then
! use loop unrolling to speed up the copying
        k = mod(length,7)
        do j = 1,k
          i = map(j)
          read_array(i) = read_array(i) + keep_array(loc-1+j)
        end do
        do j = k+1,length, 7
          i = map(j)
          read_array(i) = read_array(i) + keep_array(loc-1+j)
          i = map(j+1)
          read_array(i) = read_array(i) + keep_array(loc+j)
          i = map(j+2)
          read_array(i) = read_array(i) + keep_array(loc+j+1)
          i = map(j+3)
          read_array(i) = read_array(i) + keep_array(loc+j+2)
          i = map(j+4)
          read_array(i) = read_array(i) + keep_array(loc+j+3)
          i = map(j+5)
          read_array(i) = read_array(i) + keep_array(loc+j+4)
          i = map(j+6)
          read_array(i) = read_array(i) + keep_array(loc+j+5)
        end do

!     else if (present(add)) then
!       call daxpy(length,one,keep_array(loc),1,read_array,1)
     else
       call dcopy(length,keep_array(loc),1,read_array,1)
     end if
   end if

  end subroutine MA77_read_real
!************************************************************************

  subroutine MA77_read_discard_real(rfile,keep_array,loc,length,read_array, &
             flag,data,lp,map)

! If files are in use (rfile >= 0), call of01_read to read into read_array
! with discard = .true.
! otherwise, copy from keep_array into read_array.

   integer(short), intent(in) :: rfile ! unit for reading from
   real(wp), allocatable :: keep_array(:) ! if files not in use,
!       we copy from this array into read_array
   integer(long), intent(in) :: loc ! position in keep_array to read from
   integer(short), intent(in) :: length ! number of reals to be read
   real(wp), dimension(:), intent(inout) :: read_array ! read into this array
   integer(short), intent(out) :: flag ! error flag
   type (of01_rdata), intent(inout) :: data
   integer(short), intent(in) :: lp ! unit for of01 error messages
   integer(short), optional, intent (in) :: map(length)

! Local scalars
   integer(short) :: i,j,k
   logical :: discard  ! of01 discard parameter. Set to .true.

   flag = 0
   if (length == 0) return

   discard = .true.

   if (rfile >= 0) then
! read using of01
     if (present(map)) then
       call of01_read(rfile,loc,length,read_array,flag,data,  &
            lp=lp,map=map,discard=discard)
     else
       call of01_read(rfile,loc,length,read_array,flag,data, &
            lp=lp,discard=discard)
     end if

   else
! Copy from keep_array into read_array
     if (present(map)) then
! use loop unrolling to speed up the copying
        k = mod(length,7)
        do j = 1,k
          i = map(j)
          read_array(i) = read_array(i) + keep_array(loc-1+j)
        end do
        do j = k+1,length, 7
          i = map(j)
          read_array(i) = read_array(i) + keep_array(loc-1+j)
          i = map(j+1)
          read_array(i) = read_array(i) + keep_array(loc+j)
          i = map(j+2)
          read_array(i) = read_array(i) + keep_array(loc+j+1)
          i = map(j+3)
          read_array(i) = read_array(i) + keep_array(loc+j+2)
          i = map(j+4)
          read_array(i) = read_array(i) + keep_array(loc+j+3)
          i = map(j+5)
          read_array(i) = read_array(i) + keep_array(loc+j+4)
          i = map(j+6)
          read_array(i) = read_array(i) + keep_array(loc+j+5)
        end do

     else
       call dcopy(length,keep_array(loc),1,read_array,1)
     end if
   end if

  end subroutine MA77_read_discard_real

!************************************************************************

  subroutine MA77_write_real(rfile,size_array,keep_array,loc,length, &
             write_array,flag,data,lp,maxstore,used,inactive)

! If files are in use (rfile > 0), call of01_write to write from write_array.
! otherwise, try and copy from write_array to keep_array. If
! keep_array is not long enough, try and reallocate. Otherwise,
! switch to using a file.

   integer(short), intent(inout) :: rfile
!          index for writing to file or < 0 if in core
   real(wp), allocatable :: keep_array(:)
!         if files not in use, we copy from write_array into this array.
   integer(long), intent(inout) :: size_array ! size of keep_array
   integer(long), intent(in) :: loc     ! position in keep_array to write to
   integer(short), intent(in) :: length ! number of reals to be written
   real(wp), intent(in) :: write_array(:) ! write from this array
   integer(short), intent(out) :: flag    ! error flag
   type (of01_rdata), intent(inout) :: data
   integer(short), intent(in) :: lp       ! unit for of01 error messages
   integer(long), intent(in) :: maxstore ! max storage available
   integer(long), intent(inout) :: used  ! storage used
   integer(long), optional, intent(in) :: inactive  ! of01 inactive parameter

   real(wp), allocatable :: temp(:) ! temporary array

! Local scalars
   integer(short):: i,lw1,nwrite,st
   integer(long) :: klong,lkeep,loc1,lw,lw1_long

   flag = 0
   if (length == 0) return

   if (rfile > 0) then
! write using of01
     if (present(inactive)) then
       call of01_write(rfile,loc,length,write_array,flag,data, &
            lp=lp,inactive=inactive)
     else
       call of01_write(rfile,loc,length,write_array,flag,data,lp=lp)
     end if
     return
   end if

! Try and copy from write_array into keep_array
   lkeep = size_array
   if (loc+length-1 <= lkeep) then
     call dcopy(length,write_array,1,keep_array(loc),1)
     return
   end if

! keep_array is not long enough ... try and reallocate provided the
! limit maxstore is not exceeded.
   lw = loc - 1
   klong = max(2*lkeep,lw+length)
   if ((used+(klong-lkeep)*2) < maxstore) then
! there is space to deallocate keep_array and reallocate with length klong
     allocate (temp(lw),stat=st)
     if (st == 0) then
! take a temporary copy of the contents of keep_array and then
! reallocate keep_array with longer length
! Remember that lw is a long integer so split this copy into chunks
       if (lw > 0) then
         lw1 = lup1
         nwrite = lw/lw1
         loc1 = 1
         do i = 1,nwrite
           call dcopy(lw1,keep_array(loc1),1,temp(loc1),1)
           loc1 = loc1 + lw1
         end do
         lw1_long = lw1
         lw1 = int((lw - nwrite*lw1_long),short)
         call dcopy(lw1,keep_array(loc1),1,temp(loc1),1)
       end if

       deallocate (keep_array,stat=st)
       allocate (keep_array(klong),stat=st)
       if (st == 0) then
         size_array = klong
         if (lw > 0) then
           lw1 = lup1
           loc1 = 1
           do i = 1,nwrite
             call dcopy(lw1,temp(loc1),1,keep_array(loc1),1)
             loc1 = loc1 + lw1
           end do
           lw1_long = lw1
           lw1 = int((lw - nwrite*lw1_long),short)
           call dcopy(lw1,temp(loc1),1,keep_array(loc1),1)
         end if

         call dcopy(length,write_array,1,keep_array(loc),1)
         deallocate (temp,stat=st)
         used = used + (klong - lkeep)*2
         return
       else
! contents of keep_array are in temp but not enough room
! to allocate larger keep_array and copy back.
! Switch to using file. First write out contents of temp.
!        write (6,*) 'switch to rfile 1', rfile
         rfile = -rfile
         if (lw > 0) then
           lw1 = lup1
           nwrite = lw/lw1
           loc1 = 1
           do i = 1,nwrite
             call of01_write(rfile,loc1,lw1,temp(loc1:loc1+lw1-1), &
                  flag,data,lp=lp)
             if (flag < 0) then
               deallocate (temp,stat=st);return
             end if
             loc1 = loc1 + lw1
           end do
           lw1_long = lw1
           lw1 = int((lw - nwrite*lw1_long),short)
           call of01_write(rfile,loc1,lw1,temp(loc1:loc1+lw1-1),  &
                flag,data,lp=lp)
         end if

         if (flag >= 0) then
           if (present(inactive)) then
             call of01_write(rfile,loc,length,write_array,flag,data, &
                 lp=lp,inactive=inactive)
           else
             call of01_write(rfile,loc,length,write_array,flag,data,lp=lp)
           end if
         end if
! freed up the space used by keep_array
         used = used - lkeep*2
         deallocate (temp,stat=st)
         return
       end if
     end if
     deallocate (temp,stat=st)
   end if

! Insufficient space to write to keep_array so we must switch to
! using a file with index -rfile (we have already opened a file
! with this index). First write out contents of keep_array
! and then write_array.
!   write (6,*) 'switch to rfile 2', rfile
   rfile = -rfile
   if (lw > 0) then
     loc1 = 1
     lw1 = lup1
     nwrite = lw/lw1
     do i = 1,nwrite
       call of01_write(rfile,loc1,lw1,keep_array(loc1:loc1+lw1-1), &
            flag,data,lp=lp)
       if (flag < 0) return
       loc1 = loc1 + lw1
     end do
     lw1_long = lw1
     lw1 = int((lw - nwrite*lw1_long),short)
     call of01_write(rfile,loc1,lw1,keep_array(loc1:loc1+lw1-1), &
          flag,data,lp=lp)
     if (flag < 0) return
   end if

   if (present(inactive)) then
     call of01_write(rfile,loc,length,write_array,flag,data, &
          lp=lp,inactive=inactive)
   else
     call of01_write(rfile,loc,length,write_array,flag,data,lp=lp)
   end if

   deallocate (keep_array,stat=st)
   size_array = 0
! freed up the space used by keep_array
   used = used - lkeep*2

  end subroutine MA77_write_real

!************************************************************************

  subroutine MA77_write_integer(ifile,size_array,keep_array,loc,length, &
             write_array,flag,data,lp,maxstore,used)

! If files are in use (ifile > 0), call of01_write to write from write_array.
! otherwise, try and copy from write_array to keep_array. If
! keep_array is not long enough, try and reallocate. Otherwise,
! switch to using a file on unit unit.

   integer(short), intent(inout) :: ifile
!          index for writing to file or < 0 if in core
   integer(short), allocatable :: keep_array(:)
!          if files not in use, we copy from write_array into this array.
   integer(long), intent(inout) :: size_array ! size of keep_array
   integer(long), intent(in) :: loc ! position in keep_array to write to
   integer(short), intent(in) :: length   ! number of reals to be written
   integer(short), intent(in) :: write_array(:) ! write from this array
   integer(short), intent(out) :: flag    ! error flag
   type (of01_idata), intent(inout) :: data
   integer(short), intent(in) :: lp       ! unit for of01 error messages
   integer(long), intent(in) :: maxstore  ! max storage available
   integer(long), intent(inout) :: used   ! storage used

   integer(short), allocatable :: temp(:) ! temporary array

! Local scalars
   integer(short) :: i,k,lw1,nwrite,st
   integer(long) :: klong,kk,lkeep,lw1_long,loc1,lw

   flag = 0
   if (length == 0) return

   if (ifile > 0) then
!   write (6,*) 'ifile,loc,length,loc+length-1',&
!   ifile,loc,length,loc+length-1
! write using of01
     call of01_write(ifile,loc,length,write_array,flag,data,lp=lp)
     return
   end if

! Try and copy from write_array into keep_array
   lkeep = size_array
!  write (6,*) 'ifile,loc,length,loc+length-1,lkeep',&
!  ifile,loc,length,loc+length-1,size(keep_array)
   if (loc+length-1 <= lkeep) then
     k = mod(length,7)
     keep_array(loc:loc-1+k) = write_array(1:k)
     do i = k+1,length, 7
       keep_array(loc-1+i) = write_array(i)
       keep_array(loc+i)   = write_array(i+1)
       keep_array(loc+i+1) = write_array(i+2)
       keep_array(loc+i+2) = write_array(i+3)
       keep_array(loc+i+3) = write_array(i+4)
       keep_array(loc+i+4) = write_array(i+5)
       keep_array(loc+i+5) = write_array(i+6)
     end do
!    keep_array(loc:loc+length-1) = write_array(1:length)
     return
   end if

! keep_array is not long enough ... try and reallocate with
! space for an extra lkeep integers (provided the
! limit maxstore is not exceeded).
!  write (6,*) 'reallocate. ifile',ifile,loc+length-1,lkeep
   lw = loc - 1
   kk = max(2*lkeep,lw+length)
!  write (6,*) 'lw,used,lkeep',lw,used,lkeep
   if ((used+(kk-lkeep)) < maxstore) then
     allocate (temp(lw),stat=st)
     if (st == 0) then
! take a temporary copy of the contents of keep_array and then
! reallocate keep_array with longer length
       klong = mod(lw,7_long)
       temp(1:klong) = keep_array(1:klong)
       do i = klong+1,lw, 7
         temp(i) = keep_array(i)
         temp(i+1) = keep_array(i+1)
         temp(i+2) = keep_array(i+2)
         temp(i+3) = keep_array(i+3)
         temp(i+4) = keep_array(i+4)
         temp(i+5) = keep_array(i+5)
         temp(i+6) = keep_array(i+6)
       end do
!      temp(1:lw) = keep_array(1:lw)
       deallocate (keep_array,stat=st)

       allocate (keep_array(kk),stat=st)
       if (st == 0) then
         size_array = kk
         klong = mod(lw,7_long)
         keep_array(1:klong) = temp(1:klong)
         do i = klong+1,lw, 7_long
           keep_array(i) = temp(i)
           keep_array(i+1) = temp(i+1)
           keep_array(i+2) = temp(i+2)
           keep_array(i+3) = temp(i+3)
           keep_array(i+4) = temp(i+4)
           keep_array(i+5) = temp(i+5)
           keep_array(i+6) = temp(i+6)
         end do
!        keep_array(1:lw) = temp(1:lw)

         k = mod(length,7)
         keep_array(loc:loc-1+k) = write_array(1:k)
         do i = k+1,length, 7
           keep_array(loc-1+i) = write_array(i)
           keep_array(loc+i)   = write_array(i+1)
           keep_array(loc+i+1) = write_array(i+2)
           keep_array(loc+i+2) = write_array(i+3)
           keep_array(loc+i+3) = write_array(i+4)
           keep_array(loc+i+4) = write_array(i+5)
           keep_array(loc+i+5) = write_array(i+6)
         end do
!        keep_array(loc:loc+length-1) = write_array(1:length)

         deallocate (temp,stat=st)
         used = used + (kk - lkeep)
         return
       else
! contents of keep_array are in temp but not enough room
! to allocate larger keep_array and copy back.
! Switch to using file. First write out contents of temp.
!        write (6,*) 'switch to ifile 1', ifile
         ifile = -ifile
         if (lw > 0) then
           loc1 = 1
           lw1 = lup1
           nwrite = lw/lw1
           do i = 1,nwrite
             call of01_write(ifile,loc1,lw1,temp(loc1:loc1+lw1-1), &
                  flag,data,lp=lp)
             if (flag < 0) then
               deallocate (temp,stat=st); return
             end if
             loc1 = loc1 + lw1
           end do
           lw1_long = lw1
           lw1 = int((lw - nwrite*lw1_long),short)  
           call of01_write(ifile,loc1,lw1,temp(loc1:loc1+lw1-1), &
                flag,data,lp=lp)
         end if

         if (flag >= 0) &
           call of01_write(ifile,loc,length,write_array,flag,data,lp=lp)
         used = used - lkeep
         deallocate (temp,stat=st)
         return
       end if
     end if
     deallocate (temp,stat=st)
   end if

! Insufficient space to write to keep_array so we must switch to
! using a file with index -ifile (we have already opened a file
! with this index). First write out contents of keep_array(1:lw)
! and then write_array.
! Note: lw is long integer so need to write in short integer chunks
!   write (6,*) 'switch to ifile 2', ifile
   ifile = -ifile
   if (lw > 0) then
     loc1 = 1
     lw1 = lup1
     nwrite = lw/lw1
     do i = 1,nwrite
       call of01_write(ifile,loc1,lw1,keep_array(loc1:loc1+lw1-1), &
            flag,data,lp=lp)
       if (flag < 0) return
       loc1 = loc1 + lw1
     end do
     lw1_long = lw1
     lw1 = int((lw - nwrite*lw1_long),short)
     call of01_write(ifile,loc1,lw1,keep_array(loc1:loc1+lw1-1), &
          flag,data,lp=lp)
     if (flag < 0) return
   end if

   call of01_write(ifile,loc,length,write_array,flag,data,lp=lp)
   deallocate (keep_array,stat=st)
   size_array = 0
   used = used - lkeep

  end subroutine MA77_write_integer

!************************************************************************

  subroutine MA77_finalise_double(keep,control,info,restart_file)

! This routine must be called after all other calls to routines
! in the package.
! Calls are made to of01_close and then to of01_end.
! If restart_file is not present, components of keep are deallocated
! and keep%status reset to 0 on successful return.
! Otherwise, data is written to a sequential access file
! so that future further solves are possible.

    type (MA77_keep), intent (inout) :: keep    ! See derived-type declaration
    type (MA77_control), intent (in) :: control ! See derived-type declaration
    type (MA77_info), intent (out) :: info      ! See derived-type declaration

    character (len=*), optional, intent (in) :: restart_file
! If present, must hold the name of the file to which data is written so
! that future solves are possible.

     integer(short), allocatable :: nchild(:) ! used to hold the
!             number of children at each node of the tree.
     integer(short), allocatable :: nelim(:) ! used to hold the

    integer(long) :: lenw ! temporary variable.
    integer(long) :: ilenw ! length in pages of the part of
! the main integer virtual array that have been written.
    integer(long) :: rlenw ! length in pages of the part of
! the main reals virtual array that have been written.
    integer(short) :: flag ! error flag
    logical :: ex ! exist variable for inquire
    integer(short) :: i
    integer(short) :: ios
    integer(short) :: nout
    logical :: open
    integer(short) :: st     ! stat parameter
    integer(short) :: unit   ! unit number for file
    logical :: lkeep  ! Set to .true. if files are to be kept on closure

! Possible error returns:
!  -1   Allocation error
!  -3   Call with restart_file present but factorisation not performed
!  -5   Error from of01_close (error in Fortran inquire)
!  -7   Error from of01_close (error in Fortran open)
!  -8   Deallocation error
! -12   A file with name restart_file already exists
! -13   len(restart_file) > 500
! -15   Error in Fortran write
! -41   Error in Fortran CLOSE

! immediate return if status = 0 (nothing to be done)
!   write (6,*) 'final',keep%status
!   if (present(restart_file)) write (6,*) 'restart_file  ',restart_file
    if (keep%status == 0) return

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_finalise'
    if (control%print_level >= 2 .and. control%unit_diagnostics >= 0) then
      if (present(restart_file)) then
       write (control%unit_diagnostics,'(/a)') &
     ' Entering MA77_finalise, with restart_file present'
      else
       write (control%unit_diagnostics,'(/a)') &
     ' Entering MA77_finalise'
      end if
     end if

    st = 0

! If error occurred in MA77_open before of01 was called, all
! we have to do is deallocate components of keep
    if (keep%status == -2) then
      call final_deallocate
      keep%status = 0
      return
    end if

    lkeep = .false.
    if (present(restart_file)) then
! Perform quick checks. If an error occurs at this point, we would
! want the user to be able to recall MA77_finalise so we do no deallocation
! before returning.
! Check factorisation was performed.
      if (keep%status /= 3) then
        info%flag = -3
        call MA77_print_iflag(keep,nout,info%flag)
        if (nout >= 0) write (nout,'(a)') &
      ' restart_file is present but factors not computed'
        return
      end if
! Check restart_file.
      if (len(restart_file) > 500) then
        info%flag = -13
        call ma77_print_iflag(keep,nout,info%flag)
        return
      end if
      inquire (file=trim(restart_file),exist=ex)
      if (ex) then
        info%flag = -12
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if
! find suitable unit number on which to open file
      do unit = 10, huge(0)
        inquire (unit=unit,iostat=ios,opened=open)
        if (ios /= 0) then
! Error in inquire
          info%flag = -5; info%iostat = ios
          call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
          return
        end if
        if (.not. open) exit
      end do
! Open an unformatted sequential access file
      info%unit_restart = unit
      open (unit=unit,file=trim(restart_file),iostat=ios,status='new', &
            form='unformatted')
      if (ios /= 0 ) then
        info%iostat = ios
        info%flag = -7
        call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
        return
      end if
! main integer and real files to be kept on closure
      lkeep = .true.
    end if

! Write what is left in the buffers to the appropriate
! direct access files and close the files
    call of01_close(abs(keep%index(1)),ilenw,info%num_file(1),flag, &
         keep%idata,lp=-1,lkeep=lkeep)
    if (flag < 0) then
      info%flag = flag
      if (flag == -14) info%flag = -41
      info%iostat = keep%idata%iostat
    end if

    call of01_close(abs(keep%index(2)),rlenw,info%num_file(2),flag, &
         keep%rdata,lp=-1,lkeep=lkeep)
    if (flag < 0) then
      info%flag = flag
      if (flag == -14) info%flag = -41
      info%iostat = keep%rdata%iostat
    end if

! We do not want to keep the workspace files so close and delete.
    call of01_close(abs(keep%index(3)),lenw,info%num_file(3),flag, &
         keep%rdata,lp=-1,lkeep=.false.)
    if (flag < 0) then
     info%flag = flag
     if (flag == -14) info%flag = -41
     info%iostat = keep%rdata%iostat
    end if

    call of01_close(abs(keep%index(4)),lenw,info%num_file(4),flag, &
         keep%rdata,lp=-1,lkeep=.false.)
    if (flag < 0) then
      info%flag = flag
      if (flag == -14) info%flag = -41
      info%iostat = keep%rdata%iostat
    end if

!    write (*,*) 'info%num_file ', info%num_file(1:4)

! check whether an error was encountered in any of the above closure calls.
    if (info%flag < 0) then
!      write (6,*) 'info%flag =',info%flag
      call error_finalise
      return
    end if

! If restart is requested,
! write the integers in keep that we have to preserve to the file.
    if (present(restart_file)) then
! First preserve data that is needed by of01
      write (unit=unit,iostat=ios,err=10) &
        keep%lpage,       &
        keep%npage,       &
        keep%file_size,   &
        ilenw,            &
        rlenw,            &
        keep%index(1),    &
        keep%index(2),    &
        keep%n

 10   if (ios /= 0) then
        info%flag = -15; info%iostat = ios
        call error_finalise
        return
      end if

      if (keep%n > 0) then
        call inner_write(keep,keep%tree)
        if (info%flag < 0) then
          call error_finalise
          return
        end if
      end if
! finished writing everything we need to restart. Close the file.
      close (unit=unit,status='keep')
    end if

! deallocate components of the of01 derived types
    call of01_end(flag,keep%idata,-1)
    if (flag < 0) then
      info%flag = flag
      info%stat = keep%idata%stat
    end if
    call of01_end(flag,keep%rdata,-1)
    if (flag < 0) then
      info%flag = flag
      info%stat = keep%rdata%stat
    end if
    if (flag < 0) then
      info%flag = flag
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
    end if

! deallocate the components of tree that are finished with
     if (keep%status > 1 .and. keep%n > 0) then
       do i = keep%l1,keep%l2
         deallocate (keep%tree(i)%child,stat=st)
      end do
    end if
    deallocate (keep%tree,stat=st)

    call final_deallocate

! Set status parameter
    keep%status = 0
    if (info%flag < 0) keep%status = -1

  contains
!*******************************************
    subroutine error_finalise
      integer(short) :: flag

      call MA77_print_iflag(keep,nout,info%flag)
      keep%status = -1

      call of01_end(flag,keep%idata,-1)
      info%stat = keep%idata%stat

      call of01_end(flag,keep%rdata,-1)
      info%stat = keep%rdata%stat

! deallocate all array components of keep
      deallocate (keep%arow,stat=st)
      deallocate (keep%aelt,stat=st)
      deallocate (keep%clist,stat=st)
      deallocate (keep%iptr,stat=st)
      deallocate (keep%ifile,stat=st)
      deallocate (keep%map,stat=st)
      deallocate (keep%new,stat=st)
      deallocate (keep%rfile,stat=st)
      deallocate (keep%roots,stat=st)
      deallocate (keep%size,stat=st)
      deallocate (keep%size_ind,stat=st)
      deallocate (keep%svar,stat=st)
      deallocate (keep%splitp,stat=st)
      deallocate (keep%vars,stat=st)
      deallocate (keep%varflag,stat=st)

      deallocate (keep%imain,stat=st)
      deallocate (keep%rmain,stat=st)
      deallocate (keep%rwork,stat=st)
      deallocate (keep%rwdelay,stat=st)

      if (keep%status > 1) then
        do i = keep%l1,keep%l2
          deallocate (keep%tree(i)%child,stat=st)
        end do
      end if
      deallocate (keep%tree,stat=st)

    end subroutine error_finalise
!*******************************************
    subroutine final_deallocate

! deallocate all array components of keep
      deallocate (keep%arow,stat=st)
      deallocate (keep%aelt,stat=st)
      deallocate (keep%clist,stat=st)
      deallocate (keep%iptr,stat=st)
      deallocate (keep%ifile,stat=st)
      deallocate (keep%map,stat=st)
      deallocate (keep%new,stat=st)
      deallocate (keep%rfile,stat=st)
      deallocate (keep%roots,stat=st)
      deallocate (keep%size,stat=st)
      deallocate (keep%size_ind,stat=st)
      deallocate (keep%splitp,stat=st)
      deallocate (keep%svar,stat=st)
      deallocate (keep%vars,stat=st)
      deallocate (keep%varflag,stat=st)

      deallocate (keep%imain,stat=st)
      deallocate (keep%rmain,stat=st)
      deallocate (keep%rwork,stat=st)
      deallocate (keep%rwdelay,stat=st)

    end subroutine final_deallocate
!*******************************************
   subroutine inner_write(keep,tree)

    type (MA77_keep), intent (inout) :: keep
    type (MA77_node), intent (inout) :: tree(keep%l1:keep%l2)

    integer(short) :: i
    integer(short) :: ios
    integer(short) :: nroot
    integer(short) :: size_rfile
    integer(short) :: st

! Have to preserve keep%rfile(1:nelt) for possible residual calculations
      nroot = size(keep%roots)
      size_rfile = size(keep%rfile)

      keep%l1 = keep%nelt+1
      keep%l2 = keep%tnode

      write (unit=unit,iostat=ios,err=10) &
        nroot,                 &
        size_rfile,            &
        keep%l1,               &
        keep%l2,               &
        keep%ltree,            &
        keep%tnode,            &
        keep%nelt

! store the number of children that each (non-leaf) node in the tree has.
      deallocate (nchild,nelim,stat=st)
      allocate (nchild(keep%l1:keep%l2),nelim(keep%l1:keep%l2),stat=st)
      if (st /= 0) then
        info%flag = -1; return
      end if
      nchild = 0
      nelim = 0
      do i = keep%l1,keep%l2
        if (allocated(keep%tree(i)%child)) &
          nchild(i) = size(keep%tree(i)%child)
        nelim(i) = keep%tree(i)%nelim
      end do

      write (unit=unit,iostat=ios,err=10) &
        keep%flag,             &
        keep%lup,              &
        keep%nb,               &
        keep%element_input,    &
        keep%pos_def,          &
        keep%dfree,            &
        keep%ifree,            &
        keep%inelrs,           &
        keep%mx_ifree,         &
        keep%mvar,             &
        keep%mmvar,            &
        keep%maxdepth,         &
        keep%maxelim_actual,   &
        keep%maxfa,            &
        keep%maxfront,         &
        keep%maxfrontb,        &
        keep%maxlen,           &
        keep%maxstore,         &
        keep%name,             &
        keep%ntwo,             &
        keep%null,             &
        keep%posfac,           &
        keep%posint,           &
        keep%rfree,            &
        keep%rtopmx,           &
        keep%rtopdmx,          &
        keep%rfile,            &
        keep%roots,            &
        keep%scale,            &
        keep%used,             &
        keep%size(1:keep%tnode),     &
        keep%ifile(1:keep%tnode),    &
        nchild(keep%l1:keep%l2),     &
        nelim(keep%l1:keep%l2),      &
        keep%splitp(keep%l1:keep%l2)

      do i = keep%l1,keep%l2
        if (nchild(i) > 0) write (unit=unit,iostat=ios,err=10) tree(i)%child
      end do

      if (.not.keep%pos_def) then
! reuse nelim to hold keep%size_ind(keep%l1:keep%l2)
! (if we do not fo this, ifort gives a warning message about taking
!  temporay copy)
         do i = keep%l1,keep%l2
           nelim(i) =  keep%size_ind(i)
         end do
         write (unit=unit,iostat=ios,err=10) nelim(keep%l1:keep%l2)
      end if
      deallocate (nchild,stat=st)
      deallocate (nelim,stat=st)

      if (keep%index(1) < 0) &
! Integers are in imain
        write (unit=unit,iostat=ios,err=10) keep%imain(1:keep%ifree-1)
      if (keep%index(2) < 0) &
! Reals are in rmain
        write (unit=unit,iostat=ios,err=10) keep%rmain(1:keep%rfree-1)

 10   if (ios /= 0) then
        info%flag = -15; info%iostat = ios; return
      end if


   end subroutine inner_write


  end subroutine MA77_finalise_double
!***************************************************************

  subroutine MA77_restart_double(restart_file,filename,keep,control,info,path)

! This routine must be called if the user wants to restart, that is, perform
! further solves after a call to MA77_finalise.
! Data is written to the sequential access file restart_file.
! If direct access files were used for the factors,
! calls are made to of01_initialize and then to of01_open.

    type (MA77_keep), intent (out) :: keep      ! See derived-type declaration
    type (MA77_control), intent (in) :: control ! See derived-type declaration
    type (MA77_info), intent (inout) :: info    ! See derived-type declaration

    character (len=*), intent (in) :: restart_file
! Must hold the name of the file to which data was written by MA77_finalise.
    character (len=*), optional, intent (in) :: path(:)
    character (len=*), intent (in) :: filename(4)
! Must hold the name of the paths/files

     integer(short), allocatable :: nchild(:) ! used to hold the number of
!        children at each node of the tree.
     integer(short), allocatable :: nelim(:) ! used to hold the number of

    integer(long) :: ilenw ! length in pages of the part of
! the main integer virtual array that have been written.
    integer(long) :: rlenw ! length in pages of the part of
! the main real virtual array that have been written.
    logical :: ex ! exist variable for inquire
    integer(short) :: file
    integer(short) :: ind ! of01 index
    integer(short) :: ios
    integer(short) :: nout
    integer(short) :: nroot ! size of keep%roots
    logical :: open
    integer(short) :: size_rfile ! size of keep%rfile
    integer(short) :: st     ! stat parameter
    integer(short) :: unit   ! unit number for file

! Possible error returns:
!  -1   Allocation error
!  -5   Error in Fortran inquire
!  -6   Error in Fortran read
!  -7   Error in Fortran open
!  -8   Deallocation error
! -12   File with name filename(3:5) already exists
! -13   len(filename) > 400
! -16   len(path) > 400 or open unsuccessful for all elements of path
! -28   Files with names restart_file and/or filename(1:2) do not exist.

    nout = control%unit_error
    if (control%print_level < 0) nout = -1
    keep%name = 'MA77_restart'

    st = 0
    info%flag = 0

    if (present(path)) then
      if (len(path) > 400) then
        info%flag = -16
        call MA77_print_iflag(keep,nout,info%flag)
        return
      end if
    end if
    if (len(filename) > 400) then
      info%flag = -13
      call MA77_print_iflag(keep,nout,info%flag)
      return
    end if

! Check restart file
    inquire (file=trim(restart_file),exist=ex)
    if (.not. ex) then
      info%flag = -28
      call MA77_print_iflag(keep,nout,info%flag)
!      if (nout >= 0) write (nout,*) ' restart_file = ',trim(restart_file)
      return
    end if
! find suitable unit number on which to open file
      do unit = 10, huge(0)
        inquire (unit=unit,iostat=ios,opened=open)
        if (ios /= 0) then
! Error in inquire
          info%flag = -5; info%iostat = ios
          call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
          return
        end if
        if (.not. open) exit
      end do
      info%unit_restart = unit
! Open the sequential access file restart_file
      open (unit=unit,file=trim(restart_file),iostat=ios,status='old', &
            form='unformatted')
      if (ios /= 0) then
        info%iostat = ios
        info%flag = -7
        call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
        close (unit=unit,status='keep')
        return
      end if

! Read in of01 data
    read (unit=unit,iostat=ios,err=10) &
        keep%lpage,            &
        keep%npage,            &
        keep%file_size,        &
        ilenw,                 &
        rlenw,                 &
        keep%index(1),         &
        keep%index(2),         &
        keep%n

    if (present(path)) then
      call of01_initialize(info%flag,keep%idata,path=path,npage=keep%npage(1),&
           lpage=keep%lpage(1),file_size=keep%file_size,lp=-1)
    else
      call of01_initialize(info%flag,keep%idata,npage=keep%npage(1), &
           lpage=keep%lpage(1),file_size=keep%file_size,lp=-1)
    end if
    if (info%flag == -1) then
      info%stat = keep%idata%stat
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      close (unit=unit,status='keep')
      return
    end if

    if (present(path)) then
      call of01_initialize(info%flag,keep%rdata,path=path,npage=keep%npage(2),&
           lpage=keep%lpage(2),file_size=keep%file_size,lp=-1)
    else
      call of01_initialize(info%flag,keep%rdata,npage=keep%npage(2), &
           lpage=keep%lpage(2),file_size=keep%file_size,lp=-1)
    end if
    if (info%flag == -1) then
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      info%stat = keep%rdata%stat
      close (unit=unit,status='keep')
      return
    end if

! Open direct access files with of01.
! If files were used for the factorization, the first two contain saved
! integer and real data, respectively.
    if (keep%index(1) > 0 .and. keep%n > 0) then
      call of01_open(filename(1),ind,info%flag,keep%idata,lenw=ilenw,lp=-1)
      if (info%flag >= 0) keep%index(1) = ind
    else
      call of01_open(filename(1),ind,info%flag,keep%idata,lp=-1)
      if (info%flag >= 0)keep%index(1) = -ind
    end if
    if (info%flag < 0) then
      info%iostat = keep%idata%iostat; file = 1
      call restart_error;  return
    end if

    if (keep%index(2) > 0 .and. keep%n > 0) then
      call of01_open(filename(2),ind,info%flag,keep%rdata,lenw=rlenw,lp=-1)
      if (info%flag >= 0)keep%index(2) = ind
    else
      call of01_open(filename(2),ind,info%flag,keep%rdata,lp=-1)
      if (info%flag >= 0)keep%index(2) = -ind
    end if
    if (info%flag < 0) then
      info%iostat = keep%rdata%iostat; file = 2
      call restart_error;  return
    end if

    call of01_open(filename(3),keep%index(3),info%flag,keep%rdata,lp=-1)
    if (info%flag < 0) then
      info%iostat = keep%rdata%iostat; file = 3
      call restart_error;  return
    end if

! Secondary work array needed in the indefinite case
    call of01_open(filename(4),keep%index(4),info%flag,keep%rdata,lp=-1)
    if (info%flag < 0) then
      info%iostat = keep%rdata%iostat; file = 4
      call restart_error;  return
    end if

    info%index(1:4) = keep%index(1:4)
    if (keep%n == 0) go to 30

    read (unit=unit,iostat=ios,err=10) &
        nroot,                 &
        size_rfile,            &
        keep%l1,               &
        keep%l2,               &
        keep%ltree,            &
        keep%tnode,            &
        keep%nelt

    deallocate (keep%tree,stat=st)
    allocate (keep%tree(keep%l1:keep%l2),stat=st)
    if (st /= 0) go to 20

    deallocate (keep%roots,stat=st)
    allocate (keep%roots(nroot),stat=st)
    if (st /= 0) go to 20

    deallocate (keep%rfile,stat=st)
    allocate (keep%rfile(size_rfile),stat=st)
    if (st /= 0) go to 20

    deallocate (keep%size,stat=st)
    allocate (keep%size(1+keep%tnode),stat=st)
    if (st /= 0) go to 20

    deallocate (keep%ifile,stat=st)
    allocate (keep%ifile(1+keep%tnode),stat=st)
    if (st /= 0) go to 20

! allocate space for the number of children that each
!(non-leaf) node in the tree has.
      deallocate (nchild,nelim,stat=st)
      allocate (nchild(keep%l1:keep%l2),nelim(keep%l1:keep%l2),stat=st)
      if (st /= 0) go to 20

! allocate space for the splitp points.
      deallocate (keep%splitp,stat=st)
      allocate (keep%splitp(keep%l1:keep%l2),stat=st)
      if (st /= 0) go to 20

! Read keep plus read into tree
      call inner_read(keep,keep%tree)

      deallocate (nchild,stat=st)
      deallocate (nelim,stat=st)
      if (info%flag < 0) then
        keep%status = -1
        return
      end if

      if (keep%index(1) < 0) then
! Read integers into imain
        deallocate (keep%imain,stat=st)
        allocate (keep%imain(1:keep%ifree-1),stat=st)
        if (st /= 0) go to 20
        read (unit=unit,iostat=ios,err=10) keep%imain(1:keep%ifree-1)
        keep%size_imain = keep%ifree-1
      end if
      if (keep%index(2) < 0) then
! Read reals into rmain
        deallocate (keep%rmain,stat=st)
        allocate (keep%rmain(1:keep%rfree-1),stat=st)
        if (st /= 0) go to 20
        read (unit=unit,iostat=ios,err=10) keep%rmain(1:keep%rfree-1)
        keep%size_rmain = keep%rfree-1
      end if
      if (keep%index(3) < 0) then
! Allocate rwork
        deallocate (keep%rwork,stat=st)
        allocate (keep%rwork(1:keep%rtopmx),stat=st)
        if (st /= 0) go to 20
        keep%size_rwork = keep%rtopmx
      end if
      if (keep%index(4) < 0) then
! Allocate rwdelay
        deallocate (keep%rwdelay,stat=st)
        allocate (keep%rwdelay(1:keep%rtopdmx),stat=st)
        if (st /= 0) go to 20
        keep%size_rwdelay = keep%rtopdmx
      end if

 10 if (ios /= 0) then
      info%flag = -6
      info%iostat = ios
      call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
      deallocate (nchild,stat=st)
      deallocate (nelim,stat=st)
      keep%status = -1
      return
    end if

 20 if (st /= 0) then
      info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      deallocate (nchild,stat=st)
      deallocate (nelim,stat=st)
      keep%status = -1
      return
    end if

! finished reading everything we need to restart. Close the file and delete.
 30 rewind (unit=unit)
    close (unit=unit,status='delete')

! set keep%status to show restart been called successfully
    keep%status = 3
    if (control%print_level >= 1 .and. control%unit_diagnostics >= 0) &
      write (control%unit_diagnostics,'(/a,i4/a,4i4)') &
     ' Leaving MA77_restart with error flag     info%flag  = ', info%flag, &
     ' File indices                        info%index(1:4) = ',info%index(1:4)


   contains
!*******************************************
    subroutine restart_error
      integer(short) :: flag,num_file
      integer(long) :: lenw

      keep%status = -2
      if (info%flag == -11) info%flag = -28
      if (info%flag == -17) info%flag = -16
      call MA77_print_iflag(keep,nout,info%flag)
! close any files that had already been opened
! (keep them if they hold factor data)
        if (file > 1) then
          if (keep%index(1) > 0) then
            call of01_close(keep%index(1),lenw,num_file,flag, &
                 keep%idata,lp=-1,lkeep=.true.)
          else
            call of01_close(-keep%index(1),lenw,num_file,flag, &
                 keep%idata,lp=-1,lkeep=.false.)
          end if
        end if
        if (file > 2) then
          if (keep%index(2) > 0) then
            call of01_close(keep%index(2),lenw,num_file,flag, &
                 keep%rdata,lp=-1,lkeep=.true.)
          else
            call of01_close(-keep%index(2),lenw,num_file,flag, &
                 keep%rdata,lp=-1,lkeep=.false.)
          end if
        end if
        if (file > 3) call of01_close(keep%index(3),lenw,num_file,flag, &
           keep%rdata,lp=-1,lkeep=.false.)
        if (file > 4) call of01_close(keep%index(4),lenw,num_file,flag, &
           keep%idata,lp=-1,lkeep=.false.)

        close (unit=unit,status='keep')

    end subroutine restart_error

!*******************************************
    subroutine inner_read(keep,tree)

    type (MA77_keep), intent (inout) :: keep
    type (MA77_node), intent (inout) :: tree(keep%l1:keep%l2)

    integer(short) :: i
    integer(short) :: ios
    integer(short) :: nc
    integer(short) :: st

      read (unit=unit,iostat=ios,err=10) &
        keep%flag,             &
        keep%lup,              &
        keep%nb,               &
        keep%element_input,    &
        keep%pos_def,          &
        keep%dfree,            &
        keep%ifree,            &
        keep%inelrs,           &
        keep%mx_ifree,         &
        keep%mvar,             &
        keep%mmvar,            &
        keep%maxdepth,         &
        keep%maxelim_actual,   &
        keep%maxfa,            &
        keep%maxfront,         &
        keep%maxfrontb,        &
        keep%maxlen,           &
        keep%maxstore,         &
        keep%name,             &
        keep%ntwo,             &
        keep%null,             &
        keep%posfac,           &
        keep%posint,           &
        keep%rfree,            &
        keep%rtopmx,           &
        keep%rtopdmx,          &
        keep%rfile,            &
        keep%roots,            &
        keep%scale,            &
        keep%used,             &
        keep%size(1:keep%tnode),  &
        keep%ifile(1:keep%tnode), &
        nchild(keep%l1:keep%l2),  &
        nelim(keep%l1:keep%l2),   &
        keep%splitp(keep%l1:keep%l2)

      do i = keep%l1,keep%l2
        tree(i)%nelim = nelim(i)
        nc = nchild(i)
        allocate (tree(i)%child(nc),stat=st)
        if (st /= 0) go to 20
        if (nc > 0) read (unit=unit,iostat=ios,err=10) tree(i)%child
      end do

      if (.not.keep%pos_def) then
        allocate (keep%size_ind(keep%l1:keep%l2),stat=st)
        if (st /= 0) go to 20
! make temporary use of nelim to avoid warning message from ifort about
! making temporary copy
        read (unit=unit,iostat=ios,err=10) nelim(keep%l1:keep%l2)
        keep%size_ind(keep%l1:keep%l2) = nelim(keep%l1:keep%l2)
      end if

 10 if (ios /= 0) then
      info%flag = -6
      info%iostat = ios
      call MA77_print_iflag(keep,nout,info%flag,ios=info%iostat)
      return
    end if

 20 if (st /= 0) then
      info%flag = -1
      info%stat = st
      call MA77_print_iflag(keep,nout,info%flag,st=info%stat)
      return
    end if

    end subroutine inner_read

  end subroutine MA77_restart_double

end module hsl_MA77_double
