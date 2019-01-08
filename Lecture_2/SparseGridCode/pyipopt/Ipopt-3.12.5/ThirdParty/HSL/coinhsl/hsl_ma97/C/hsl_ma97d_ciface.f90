module hsl_ma97_double_ciface
   use iso_c_binding
   use hsl_ma97_double, only:                            &
      f_ma97_akeep            => ma97_akeep,             &
      f_ma97_control          => ma97_control,           &
      f_ma97_fkeep            => ma97_fkeep,             &
      f_ma97_info             => ma97_info,              &
      f_ma97_analyse          => ma97_analyse,           &
      f_ma97_analyse_coord    => ma97_analyse_coord,     &
      f_ma97_factor           => ma97_factor,            &
      f_ma97_factor_solve     => ma97_factor_solve,      &
      f_ma97_solve            => ma97_solve,             &
      f_ma97_free             => ma97_free,              &
      f_ma97_enquire_posdef   => ma97_enquire_posdef,    &
      f_ma97_enquire_indef    => ma97_enquire_indef,     &
      f_ma97_alter            => ma97_alter,             &
      f_ma97_solve_fredholm   => ma97_solve_fredholm,    &
      f_ma97_lmultiply        => ma97_lmultiply,         &
      f_ma97_sparse_fwd_solve => ma97_sparse_fwd_solve,  &
      f_ma97_get_n__          => ma97_get_n__,           &
      f_ma97_get_nz__         => ma97_get_nz__

   integer, parameter :: wp = C_DOUBLE ! pkg type
   integer, parameter :: rp = C_DOUBLE ! real type

   type, bind(C) :: ma97_control
      integer(C_INT) :: f_arrays ! true(!=0) or false(==0)
      integer(C_INT) :: action ! true(!=0) or false(==0)
      integer(C_INT) :: nemin
      real(rp) :: multiplier
      integer(C_INT) :: ordering
      integer(C_INT) :: print_level
      integer(C_INT) :: scaling
      real(rp) :: small
      real(rp) :: u
      integer(C_INT) :: unit_diagnostics
      integer(C_INT) :: unit_error
      integer(C_INT) :: unit_warning
      integer(C_LONG) :: factor_min
      integer(C_INT) :: solve_blas3
      integer(C_LONG) :: solve_min
      integer(C_INT) :: solve_mf
      real(rp) :: consist_tol
      integer(C_INT) :: ispare(5)
      real(rp) :: rspare(10)
   end type ma97_control

   type, bind(C) :: ma97_info
      integer(C_INT) :: flag
      integer(C_INT) :: flag68
      integer(C_INT) :: flag77
      integer(C_INT) :: matrix_dup
      integer(C_INT) :: matrix_rank
      integer(C_INT) :: matrix_outrange
      integer(C_INT) :: matrix_missing_diag
      integer(C_INT) :: maxdepth
      integer(C_INT) :: maxfront
      integer(C_INT) :: num_delay
      integer(C_LONG) :: num_factor
      integer(C_LONG) :: num_flops
      integer(C_INT) :: num_neg
      integer(C_INT) :: num_sup
      integer(C_INT) :: num_two
      integer(C_INT) :: ordering
      integer(C_INT) :: stat
      integer(C_INT) :: ispare(5)
      real(rp) :: rspare(10)
   end type ma97_info
contains
   subroutine copy_control_in(ccontrol, fcontrol, f_arrays)
      type(ma97_control), intent(in) :: ccontrol
      type(f_ma97_control), intent(out) :: fcontrol
      logical, intent(out) :: f_arrays

      f_arrays                   = (ccontrol%f_arrays.ne.0)
      fcontrol%action            = (ccontrol%action.ne.0)
      fcontrol%nemin             = ccontrol%nemin
      fcontrol%multiplier        = ccontrol%multiplier
      fcontrol%ordering          = ccontrol%ordering
      fcontrol%print_level       = ccontrol%print_level
      fcontrol%scaling           = ccontrol%scaling
      fcontrol%small             = ccontrol%small
      fcontrol%u                 = ccontrol%u
      fcontrol%unit_diagnostics  = ccontrol%unit_diagnostics
      fcontrol%unit_error        = ccontrol%unit_error
      fcontrol%unit_warning      = ccontrol%unit_warning
      fcontrol%factor_min        = ccontrol%factor_min
      fcontrol%solve_blas3       = (ccontrol%solve_blas3.ne.0)
      fcontrol%solve_min         = ccontrol%solve_min
      fcontrol%solve_mf          = (ccontrol%solve_mf.ne.0)
      fcontrol%consist_tol       = ccontrol%consist_tol
   end subroutine copy_control_in

   subroutine copy_info_out(finfo,cinfo)
      type(f_ma97_info), intent(in) :: finfo
      type(ma97_info), intent(out) :: cinfo

      cinfo%flag                 = finfo%flag
      cinfo%flag68               = finfo%flag68
      cinfo%flag77               = finfo%flag77
      cinfo%matrix_dup           = finfo%matrix_dup
      cinfo%matrix_rank          = finfo%matrix_rank
      cinfo%matrix_outrange      = finfo%matrix_outrange
      cinfo%matrix_missing_diag  = finfo%matrix_missing_diag
      cinfo%maxdepth             = finfo%maxdepth
      cinfo%maxfront             = finfo%maxfront
      cinfo%num_delay            = finfo%num_delay
      cinfo%num_factor           = finfo%num_factor
      cinfo%num_flops            = finfo%num_flops
      cinfo%num_neg              = finfo%num_neg
      cinfo%num_sup              = finfo%num_sup
      cinfo%num_two              = finfo%num_two
      cinfo%ordering             = finfo%ordering
      cinfo%stat                 = finfo%stat
   end subroutine copy_info_out
end module hsl_ma97_double_ciface

subroutine ma97_default_control_d(ccontrol) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(ma97_control), intent(out) :: ccontrol

   type(f_ma97_control) :: fcontrol

   ccontrol%f_arrays = 0 ! false
   if(fcontrol%action) then
      ccontrol%action = 1 ! true
   else
      ccontrol%action = 0 ! false
   endif
   ccontrol%nemin             = fcontrol%nemin
   ccontrol%multiplier        = fcontrol%multiplier
   ccontrol%ordering          = fcontrol%ordering
   ccontrol%print_level       = fcontrol%print_level
   ccontrol%scaling           = fcontrol%scaling
   ccontrol%small             = fcontrol%small
   ccontrol%u                 = fcontrol%u
   ccontrol%unit_diagnostics  = fcontrol%unit_diagnostics
   ccontrol%unit_error        = fcontrol%unit_error
   ccontrol%unit_warning      = fcontrol%unit_warning
   ccontrol%factor_min        = fcontrol%factor_min
   if(fcontrol%solve_blas3) then
      ccontrol%solve_blas3 = 1 ! true
   else
      ccontrol%solve_blas3 = 0 ! false
   endif
   ccontrol%solve_min         = fcontrol%solve_min
   if(fcontrol%solve_mf) then
      ccontrol%solve_mf = 1 ! true
   else
      ccontrol%solve_mf = 0 ! false
   endif
   ccontrol%consist_tol       = fcontrol%consist_tol
end subroutine ma97_default_control_d

subroutine ma97_analyse_d(ccheck, n, cptr, crow, cval, cakeep, ccontrol, &
      cinfo, corder) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: ccheck
   integer(C_INT), value, intent(in) :: n
   type(C_PTR), value :: cptr
   type(C_PTR), value :: crow
   type(C_PTR), value :: cval
   type(C_PTR) :: cakeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: corder

   logical :: fcheck
   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   real(wp), dimension(:), pointer :: fval
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   fcheck = (ccheck.ne.0)
   call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
   if(.not.f_arrays) then
      allocate(fptr_alloc(n+1))
      fptr_alloc(:) = fptr(:) + 1
      fptr => fptr_alloc
   endif
   call C_F_POINTER(crow, frow, shape = (/ fptr(n+1)-1 /) )
   if(.not.f_arrays) then
      allocate(frow_alloc(fptr(n+1)-1))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape = (/ fptr(n+1)-1 /) )
   else
      nullify(fval)
   endif
   if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape = (/ n /) )
   else
      nullify(forder)
   endif
   if(.not.f_arrays .and. associated(forder)) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif

   ! Allocate space to store akeep and arrange a C pointer to it
   allocate(fakeep)
   cakeep = c_loc(fakeep)

   ! Call the Fortran routine
   if(associated(forder)) then
      if(associated(fval)) then
         call f_ma97_analyse(fcheck, n, fptr, frow, fakeep, fcontrol, &
            finfo, order=forder, val=fval)
      else
         call f_ma97_analyse(fcheck, n, fptr, frow, fakeep, fcontrol, &
            finfo, order=forder)
      endif
   else
      if(associated(fval)) then
         call f_ma97_analyse(fcheck, n, fptr, frow, fakeep, fcontrol, &
            finfo, val=fval)
      else
         call f_ma97_analyse(fcheck, n, fptr, frow, fakeep, fcontrol, &
            finfo)
      endif
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)

   ! Copy order out if using C indexing
   if(.not.f_arrays .and. associated(forder)) then
      call C_F_POINTER(corder, forder, shape = (/ n /) )
      forder(:) = forder_alloc(:) - 1
   endif
end subroutine ma97_analyse_d

subroutine ma97_analyse_coord_d(n, ne, crow, ccol, cval, cakeep, ccontrol, &
      cinfo, corder) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: n
   integer(C_INT), value, intent(in) :: ne
   type(C_PTR), value :: crow
   type(C_PTR), value :: ccol
   type(C_PTR), value :: cval
   type(C_PTR) :: cakeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: corder

   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   integer(C_INT), dimension(:), pointer :: fcol
   integer, dimension(:), allocatable, target :: fcol_alloc
   integer(C_INT), dimension(:), pointer :: forder
   integer, dimension(:), allocatable, target :: forder_alloc
   real(wp), dimension(:), pointer :: fval
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(crow, frow, shape = (/ ne /) )
   if(.not.f_arrays) then
      allocate(frow_alloc(ne))
      frow_alloc(:) = frow(:) + 1
      frow => frow_alloc
   endif
   call C_F_POINTER(ccol, fcol, shape = (/ ne /) )
   if(.not.f_arrays) then
      allocate(fcol_alloc(ne))
      fcol_alloc(:) = fcol(:) + 1
      fcol => fcol_alloc
   endif
   if(C_ASSOCIATED(cval)) then
      call C_F_POINTER(cval, fval, shape = (/ ne /) )
   else
      nullify(fval)
   endif
   if(C_ASSOCIATED(corder)) then
      call C_F_POINTER(corder, forder, shape = (/ n /) )
   else
      nullify(forder)
   endif
   if(.not.f_arrays .and. associated(forder)) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif

   ! Allocate space to store akeep and arrange a C pointer to it
   allocate(fakeep)
   cakeep = c_loc(fakeep)

   ! Call the Fortran routine
   if(associated(forder)) then
      if(associated(fval)) then
         call f_ma97_analyse_coord(n, ne, frow, fcol, fakeep, fcontrol, &
            finfo, order=forder, val=fval)
      else
         call f_ma97_analyse_coord(n, ne, frow, fcol, fakeep, fcontrol, &
            finfo, order=forder)
      endif
   else
      if(associated(fval)) then
         call f_ma97_analyse_coord(n, ne, frow, fcol, fakeep, fcontrol, &
            finfo, val=fval)
      else
         call f_ma97_analyse_coord(n, ne, frow, fcol, fakeep, fcontrol, &
            finfo)
      endif
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)

   ! Copy order out if using C indexing
   if(.not.f_arrays .and. associated(forder)) then
      call C_F_POINTER(corder, forder, shape = (/ n /) )
      forder(:) = forder_alloc(:) - 1
   endif
end subroutine ma97_analyse_coord_d

subroutine ma97_factor_d(matrix_type, cptr, crow, cval, cakeep, cfkeep, &
      ccontrol, cinfo, cscale) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: matrix_type
   type(C_PTR), value :: cptr
   type(C_PTR), value :: crow
   type(C_PTR), value :: cval
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: cscale

   integer :: n, nz
   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   real(wp), dimension(:), pointer :: fval
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   real(rp), dimension(:), pointer :: fscale
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cakeep, fakeep)
   n = f_ma97_get_n__(fakeep)
   nz = f_ma97_get_nz__(fakeep)
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   if(C_ASSOCIATED(cptr)) then
      call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
      if(.not.f_arrays) then
         allocate(fptr_alloc(n+1))
         fptr_alloc(:) = fptr(:) + 1
         fptr => fptr_alloc
      endif
   else
      nullify(fptr)
   endif
   if(C_ASSOCIATED(crow)) then
      call C_F_POINTER(crow, frow, shape = (/ nz /) )
      if(.not.f_arrays) then
         allocate(frow_alloc(nz))
         frow_alloc(:) = frow(:) + 1
         frow => frow_alloc
      endif
   else
      nullify(frow)
   endif
   call C_F_POINTER(cval, fval, shape = (/ nz /) )
   if(C_ASSOCIATED(cscale)) then
      call C_F_POINTER(cscale, fscale, shape = (/ n /))
   else
      nullify(fscale)
   endif

   if(.not.C_ASSOCIATED(cfkeep)) then
      allocate(ffkeep)
      cfkeep = C_LOC(ffkeep)
   else
      call C_F_POINTER(cfkeep, ffkeep)
   endif

   ! Call the Fortran routine
   if(associated(fscale)) then
      if(associated(fptr)) then
         if(associated(frow)) then
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, ptr=fptr, row=frow, scale=fscale)
         else
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, ptr=fptr, scale=fscale)
         endif
      else
         if(associated(frow)) then
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, row=frow, scale=fscale)
         else
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, scale=fscale)
         endif
      endif
   else
      if(associated(fptr)) then
         if(associated(frow)) then
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, ptr=fptr, row=frow)
         else
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, ptr=fptr)
         endif
      else
         if(associated(frow)) then
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo, row=frow)
         else
            call f_ma97_factor(matrix_type, fval, fakeep, ffkeep, fcontrol, &
               finfo)
         endif
      endif
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_factor_d

subroutine ma97_factor_solve_d(matrix_type, cptr, crow, cval, nrhs, cx, ldx, &
      cakeep, cfkeep, ccontrol, cinfo, cscale) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: matrix_type
   type(C_PTR), value :: cptr
   type(C_PTR), value :: crow
   type(C_PTR), value :: cval
   integer(C_INT), value, intent(in) :: nrhs
   type(C_PTR), value :: cx
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: cscale

   integer :: n, nz
   integer(C_INT), dimension(:), pointer :: fptr
   integer, dimension(:), allocatable, target :: fptr_alloc
   integer(C_INT), dimension(:), pointer :: frow
   integer, dimension(:), allocatable, target :: frow_alloc
   real(wp), dimension(:), pointer :: fval
   real(wp), dimension(:,:), pointer :: fx
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   real(rp), dimension(:), pointer :: fscale
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call C_F_POINTER(cakeep, fakeep)
   n = f_ma97_get_n__(fakeep)
   nz = f_ma97_get_nz__(fakeep)
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   if(C_ASSOCIATED(cptr)) then
      call C_F_POINTER(cptr, fptr, shape = (/ n+1 /) )
      if(.not.f_arrays) then
         allocate(fptr_alloc(n+1))
         fptr_alloc(:) = fptr(:) + 1
         fptr => fptr_alloc
      endif
   else
      nullify(fptr)
   endif
   if(C_ASSOCIATED(crow)) then
      call C_F_POINTER(crow, frow, shape = (/ nz /) )
      if(.not.f_arrays) then
         allocate(frow_alloc(nz))
         frow_alloc(:) = frow(:) + 1
         frow => frow_alloc
      endif
   else
      nullify(frow)
   endif
   call C_F_POINTER(cval, fval, shape = (/ nz /) )
   call C_F_POINTER(cx, fx, shape = (/ ldx, nrhs /) )
   if(C_ASSOCIATED(cscale)) then
      call C_F_POINTER(cscale, fscale, shape = (/ n /))
   else
      nullify(fscale)
   endif

   if(.not.C_ASSOCIATED(cfkeep)) then
      allocate(ffkeep)
      cfkeep = C_LOC(ffkeep)
   else
      call C_F_POINTER(cfkeep, ffkeep)
   endif

   ! Call the Fortran routine
   if(associated(fscale)) then
      if(associated(fptr)) then
         if(associated(frow)) then
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, ptr=fptr, row=frow, scale=fscale)
         else
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, ptr=fptr, scale=fscale)
         endif
      else
         if(associated(frow)) then
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, row=frow, scale=fscale)
         else
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, scale=fscale)
         endif
      endif
   else
      if(associated(fptr)) then
         if(associated(frow)) then
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, ptr=fptr, row=frow)
         else
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, ptr=fptr)
         endif
      else
         if(associated(frow)) then
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo, row=frow)
         else
            call f_ma97_factor_solve(matrix_type, fval, nrhs, fx, ldx, fakeep, &
               ffkeep, fcontrol, finfo)
         endif
      endif
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_factor_solve_d

subroutine ma97_solve_d(job, nrhs, cx, ldx, cakeep, cfkeep, ccontrol, cinfo) &
      bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: job
   integer(C_INT), value, intent(in) :: nrhs
   type(C_PTR), value :: cx
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo

   real(wp), dimension(:,:), pointer :: fx
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cx, fx, shape = (/ ldx, nrhs /) )
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)

   ! Call the Fortran routine
   if(job.eq.0) then
      call f_ma97_solve(nrhs, fx, ldx, fakeep, ffkeep, fcontrol, finfo)
   else
      call f_ma97_solve(nrhs, fx, ldx, fakeep, ffkeep, fcontrol, finfo, job=job)
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_solve_d

subroutine ma97_free_akeep_d(cakeep) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR) :: cakeep

   type(f_ma97_akeep), pointer :: fakeep ! May be NULL

   if(C_ASSOCIATED(cakeep)) then
      call C_F_POINTER(cakeep, fakeep)
      call f_ma97_free(fakeep)
      deallocate(fakeep)
      cakeep = C_NULL_PTR
   endif
end subroutine ma97_free_akeep_d

subroutine ma97_free_fkeep_d(cfkeep) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR) :: cfkeep

   type(f_ma97_fkeep), pointer :: ffkeep ! May be NULL

   if(C_ASSOCIATED(cfkeep)) then
      call C_F_POINTER(cfkeep, ffkeep)
      call f_ma97_free(ffkeep)
      deallocate(ffkeep)
      cfkeep = C_NULL_PTR
   endif
end subroutine ma97_free_fkeep_d

subroutine ma97_finalise_d(cakeep, cfkeep) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep

   interface
      subroutine ma97_free_akeep_d(cakeep) bind(C)
         use hsl_ma97_double_ciface
         implicit none
         type(C_PTR) :: cakeep
      end subroutine ma97_free_akeep_d
      subroutine ma97_free_fkeep_d(cfkeep) bind(C)
         use hsl_ma97_double_ciface
         implicit none
         type(C_PTR) :: cfkeep
      end subroutine ma97_free_fkeep_d
   end interface

   call ma97_free_akeep_d(cakeep)
   call ma97_free_fkeep_d(cfkeep)
end subroutine ma97_finalise_d

subroutine ma97_enquire_posdef_d(cakeep, cfkeep, ccontrol, cinfo, cd) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: cd

   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   real(rp), dimension(:), pointer :: fd
   logical :: f_arrays
   integer :: n

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cakeep, fakeep)
   n = f_ma97_get_n__(fakeep)
   call C_F_POINTER(cfkeep, ffkeep)
   call C_F_POINTER(cd, fd, shape = (/ n /) )

   ! Call the Fortran routine
   call f_ma97_enquire_posdef(fakeep, ffkeep, fcontrol, finfo, fd)

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_enquire_posdef_d

subroutine ma97_enquire_indef_d(cakeep, cfkeep, ccontrol, cinfo, cpiv_order, &
      cd) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo
   type(C_PTR), value :: cpiv_order
   type(C_PTR), value :: cd

   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   real(wp), dimension(:,:), pointer :: fd
   integer(C_INT), dimension(:), pointer :: fpiv_order
   integer(C_INT), dimension(:), allocatable, target :: fpiv_order_alloc
   logical :: f_arrays
   integer :: n

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cakeep, fakeep)
   n = f_ma97_get_n__(fakeep)
   call C_F_POINTER(cfkeep, ffkeep)
   if(C_ASSOCIATED(cpiv_order)) then
      call C_F_POINTER(cpiv_order, fpiv_order, shape = (/ n /) )
      if(.not.f_arrays) then
         allocate(fpiv_order_alloc(n))
         fpiv_order_alloc(:) = fpiv_order(:) + 1
         fpiv_order => fpiv_order_alloc
      endif
   else
      nullify(fpiv_order)
   endif
   if(C_ASSOCIATED(cd)) then
      call C_F_POINTER(cd, fd, shape = (/ n, 2 /) )
   else
      nullify(fd)
   endif

   ! Call the Fortran routine
   if(associated(fpiv_order)) then
      if(associated(fd)) then
         call f_ma97_enquire_indef(fakeep, ffkeep, fcontrol, finfo, &
            piv_order=fpiv_order, d=fd)
      else
         call f_ma97_enquire_indef(fakeep, ffkeep, fcontrol, finfo, &
            piv_order=fpiv_order)
      endif
   else
      if(associated(fd)) then
         call f_ma97_enquire_indef(fakeep, ffkeep, fcontrol, finfo, d=fd)
      else
         call f_ma97_enquire_indef(fakeep, ffkeep, fcontrol, finfo)
      endif
   endif

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)

   ! Copy piv_order out if using C indexing
   if(.not.f_arrays .and. associated(fpiv_order)) then
      call C_F_POINTER(cpiv_order, fpiv_order, shape = (/ n /) )
      fpiv_order(:) = abs(fpiv_order_alloc(:)) - 1
   endif
end subroutine ma97_enquire_indef_d

subroutine ma97_alter_d(cd, cakeep, cfkeep, ccontrol, cinfo) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   type(C_PTR), value :: cd
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo

   real(wp), dimension(:,:), pointer :: fd
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays
   integer :: n

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cakeep, fakeep)
   n = f_ma97_get_n__(fakeep)
   call C_F_POINTER(cd, fd, shape = (/ n, 2 /) )
   call C_F_POINTER(cfkeep, ffkeep)

   call f_ma97_alter(fd, fakeep, ffkeep, fcontrol, finfo)

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_alter_d

subroutine ma97_solve_fredholm_d(nrhs, cflag_out, cx, ldx, cakeep, cfkeep, &
      ccontrol, cinfo) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: nrhs
   type(C_PTR), value :: cflag_out
   type(C_PTR), value :: cx
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo

   integer(C_INT), dimension(:), pointer :: fflag_out
   real(wp), dimension(:,:), pointer :: fx
   logical, dimension(:), allocatable :: fflag_out_logical
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cflag_out, fflag_out, shape = (/ nrhs /) )
   call C_F_POINTER(cx, fx, shape = (/ ldx, 2*nrhs /) )
   allocate(fflag_out_logical(nrhs))
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)

   ! Call the Fortran routine
   call f_ma97_solve_fredholm(nrhs, fflag_out_logical, fx, ldx, &
      fakeep, ffkeep, fcontrol, finfo)

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)

   ! Copy flag_out to C version
   where(fflag_out_logical(:))
      fflag_out(:) = 1
   elsewhere
      fflag_out(:) = 0
   end where
end subroutine ma97_solve_fredholm_d

subroutine ma97_lmultiply_d(ctrans, k, cx, ldx, cy, ldy, cakeep, cfkeep, &
      ccontrol, cinfo) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: ctrans
   integer(C_INT), value, intent(in) :: k
   type(C_PTR), value :: cx
   integer(C_INT), value, intent(in) :: ldx
   type(C_PTR), value :: cy
   integer(C_INT), value, intent(in) :: ldy
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo

   logical :: ftrans
   real(wp), dimension(:,:), pointer :: fx
   real(wp), dimension(:,:), pointer :: fy
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   ftrans = (ctrans.ne.0)
   call C_F_POINTER(cx, fx, shape = (/ ldx, k /) )
   call C_F_POINTER(cy, fy, shape = (/ ldx, k /) )
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)

   ! Call the Fortran routine
   call f_ma97_lmultiply(ftrans, k, fx, ldx, fy, ldy, fakeep, ffkeep, &
      fcontrol, finfo)

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)
end subroutine ma97_lmultiply_d

subroutine ma97_sparse_fwd_solve_d(nbi, cbindex, cb, corder, cnxi, &
      cxindex, cx, cakeep, cfkeep, ccontrol, cinfo) bind(C)
   use hsl_ma97_double_ciface
   implicit none

   integer(C_INT), value, intent(in) :: nbi
   type(C_PTR), value :: cbindex
   type(C_PTR), value :: cb
   type(C_PTR), value :: corder
   type(C_PTR), value :: cnxi
   type(C_PTR), value :: cxindex
   type(C_PTR), value :: cx
   type(C_PTR) :: cakeep
   type(C_PTR) :: cfkeep
   type(ma97_control), intent(in) :: ccontrol
   type(ma97_info), intent(out) :: cinfo

   integer(C_INT), dimension(:), pointer :: fbindex
   integer(C_INT), dimension(:), allocatable, target :: fbindex_alloc
   real(wp), dimension(:), pointer :: fb
   integer(C_INT), dimension(:), pointer :: forder
   integer(C_INT), dimension(:), allocatable, target :: forder_alloc
   integer(C_INT), pointer :: fnxi
   integer(C_INT), dimension(:), pointer :: fxindex
   real(wp), dimension(:), pointer :: fx
   type(f_ma97_akeep), pointer :: fakeep
   type(f_ma97_fkeep), pointer :: ffkeep
   type(f_ma97_control) :: fcontrol
   type(f_ma97_info) :: finfo
   logical :: f_arrays
   integer :: n

   logical, dimension(:), allocatable :: lflag ! FIXME inefficient

   ! Copy data in and associate pointers correctly
   call copy_control_in(ccontrol, fcontrol, f_arrays)
   call C_F_POINTER(cakeep, fakeep)
   call C_F_POINTER(cfkeep, ffkeep)
   n = f_ma97_get_n__(fakeep)
   call C_F_POINTER(cbindex, fbindex, shape = (/ nbi /) )
   if(.not.f_arrays) then
      allocate(fbindex_alloc(nbi))
      fbindex_alloc(:) = fbindex(:) + 1
      fbindex => fbindex_alloc
   endif
   call C_F_POINTER(cb, fb, shape = (/ n /) )
   call C_F_POINTER(corder, forder, shape = (/ n /) )
   if(.not.f_arrays) then
      allocate(forder_alloc(n))
      forder_alloc(:) = forder(:) + 1
      forder => forder_alloc
   endif
   call C_F_POINTER(cnxi, fnxi)
   call C_F_POINTER(cxindex, fxindex, shape = (/ n /) )
   call C_F_POINTER(cx, fx, shape = (/ n /) )

   ! FIXME inefficient
   allocate(lflag(n))
   lflag(:) = .false.

   ! Call the Fortran routine
   call f_ma97_sparse_fwd_solve(nbi, fbindex, fb, forder, lflag, fnxi, &
      fxindex, fx, fakeep, ffkeep, fcontrol, finfo)

   ! Copy information out to the C structure
   call copy_info_out(finfo, cinfo)

   ! Modify xindex for C ordering
   if(.not.f_arrays) then
      fxindex(:) = fxindex(:) - 1
   endif
end subroutine ma97_sparse_fwd_solve_d
