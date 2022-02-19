!==============================================================================!
  subroutine Bundle_Options(opts, n_opts, c_opts)
!------------------------------------------------------------------------------!
  use Const_Mod
!----------------------------------[Modules]-----------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(SL) :: opts(128)
  integer       :: n_opts
  character(QL) :: c_opts     ! catendated options
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, lc, ln
!==============================================================================!

  c_opts = ''                       ! start with an empty string
  ln = len_trim(opts(1))            ! take its length
  c_opts(1:ln) = opts(1)            ! and append first option to the end
  do i = 2, n_opts                  ! browse through remaining options ...
    lc = len_trim(c_opts)           ! ... taking the current length
    ln = len_trim(opts(i))          ! ... lenth of the new option
    c_opts(lc+1:lc+1) = ' '         ! ... adding a space
    c_opts(lc+2:lc+2+ln) = opts(i)  ! ... and the new option
  end do

  end subroutine

!==============================================================================!
  subroutine Read_Control_Petsc_Options(Flow, turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads options for PETSc solver from control file.                          !
!   I don't put it together with Read_Control_Numerical because I am afraid    !
!   it could grow a lot as we learn new PETSc options and want more from it    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Field_Mod,    only: Field_Type
  use Var_Mod,      only: Var_Type
  use Turb_Mod,     only: Turb_Type
  use Vof_Mod
  use Control_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, p, fun
  type(Var_Type),  pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: scalar(:)
  type(Var_Type),  pointer :: tq, phi
  logical                  :: found
  character(SL)            :: sstring, pstring
  character(SL)            :: opts(128)
  character(QL)            :: c_opts     ! catendated options
  integer                  :: i, n_opts, sc
  real                     :: tol
!==============================================================================!

  ! If it wasn't compiled with PETSc, don't confuse a user with this
  if(.not. PETSC_ACTIVE) return

  ! Take aliases
  Grid   => Flow % pnt_grid
  t      => Flow % t
  p      => Flow % p
  scalar => Flow % scalar
  vis    => turb % vis
  t2     => turb % t2
  fun    => Vof % fun
  call Flow % Alias_Momentum      (u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)

  !----------------------------!
  !   For momentum equations   !
  !----------------------------!
  call Control_Mod_Position_At_One_Key('PETSC_OPTIONS_FOR_MOMENTUM',  &
                                        found,                        &
                                        .false.)
  if(found) then
    call Control_Mod_Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control_Mod_Read_Char_Item_On('PREC', 'asm', pstring, .true.)
    call Control_Mod_Read_Strings_On('PREC_OPTIONS', opts, n_opts, .false.)
    if(n_opts > 0) then
      call Bundle_Options(opts, n_opts, c_opts)
    else
      c_opts = ''
    end if
    call Control_Mod_Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)
  else
    print *, '# NOTE!  PETSc options for momentum are not specified'
    print *, '#        Using the default values'
  end if

  u % solver = sstring
  v % solver = sstring
  w % solver = sstring
  u % prec = pstring
  v % prec = pstring
  w % prec = pstring
  u % prec_opts = c_opts
  v % prec_opts = c_opts
  w % prec_opts = c_opts
  u % tol = tol
  v % tol = tol
  w % tol = tol

  !---------------------------!
  !   For pressure equation   !
  !---------------------------!
  call Control_Mod_Position_At_One_Key('PETSC_OPTIONS_FOR_PRESSURE',  &
                                        found,                        &
                                        .false.)
  if(found) then
    call Control_Mod_Read_Char_Item_On('SOLVER', 'cg', sstring, .true.)
    call Control_Mod_Read_Char_Item_On('PREC', 'asm', pstring, .true.)
    call Control_Mod_Read_Strings_On('PREC_OPTIONS', opts, n_opts, .false.)
    if(n_opts > 0) then
      call Bundle_Options(opts, n_opts, c_opts)
    else
      c_opts = ''
    end if
    call Control_Mod_Read_Real_Item_On('TOLERANCE', 1.0e-5, tol, .true.)
  else
    print *, '# NOTE!  PETSc options for pressure are not specified'
    print *, '#        Using the default values'
  end if

  Flow % pp % solver = sstring
  Flow % pp % prec = pstring
  Flow % pp % prec_opts = c_opts
  Flow % pp % tol = tol

  !-------------------------!
  !   For energy equation   !
  !-------------------------!
  call Control_Mod_Position_At_One_Key('PETSC_OPTIONS_FOR_ENERGY',  &
                                        found,                      &
                                        .false.)
  if(found) then
    call Control_Mod_Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control_Mod_Read_Char_Item_On('PREC', 'asm', pstring, .true.)
    call Control_Mod_Read_Strings_On('PREC_OPTIONS', opts, n_opts, .false.)
    if(n_opts > 0) then
      call Bundle_Options(opts, n_opts, c_opts)
    else
      c_opts = ''
    end if
    call Control_Mod_Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)
  else
    print *, '# NOTE!  PETSc options for energy are not specified'
    print *, '#        Using the default values'
  end if

  t % solver = sstring
  t % prec = pstring
  t % prec_opts = c_opts
  t % tol = tol

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  call Control_Mod_Position_At_One_Key('PETSC_OPTIONS_FOR_SCALARS',  &
                                        found,                       &
                                        .false.)

  if(found) then
    call Control_Mod_Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control_Mod_Read_Char_Item_On('PREC', 'asm', pstring, .true.)
    call Control_Mod_Read_Strings_On('PREC_OPTIONS', opts, n_opts, .false.)
    if(n_opts > 0) then
      call Bundle_Options(opts, n_opts, c_opts)
      t % prec_opts = c_opts
    else
      t % prec_opts = ''
    end if
    call Control_Mod_Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)
  else
    print *, '# NOTE!  PETSc options for energy are not specified'
    print *, '#        Using the default values'
  end if

  do sc = 1, Flow % n_scalars
    phi => Flow % scalar(sc)
    phi % solver = sstring
    phi % prec = pstring
    phi % prec_opts = c_opts
    phi % tol = tol
  end do

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  call Control_Mod_Position_At_One_Key('PETSC_OPTIONS_FOR_TURBULENCE',  &
                                        found,                          &
                                        .false.)
  if(found) then
    call Control_Mod_Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control_Mod_Read_Char_Item_On('PREC', 'asm', pstring, .true.)
    call Control_Mod_Read_Strings_On('PREC_OPTIONS', opts, n_opts, .false.)
    if(n_opts > 0) then
      call Bundle_Options(opts, n_opts, c_opts)
    else
      c_opts = ''
    end if
    call Control_Mod_Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)
  else
    print *, '# NOTE!  PETSc options for energy are not specified'
    print *, '#        Using the default values'
  end if

  do i = 1, 12
    if(i .eq.  1) tq => turb % kin
    if(i .eq.  2) tq => turb % eps
    if(i .eq.  3) tq => turb % zeta
    if(i .eq.  4) tq => turb % f22
    if(i .eq.  5) tq => turb % vis
    if(i .eq.  6) tq => turb % t2
    if(i .eq.  7) tq => turb % uu
    if(i .eq.  8) tq => turb % vv
    if(i .eq.  9) tq => turb % ww
    if(i .eq. 10) tq => turb % uv
    if(i .eq. 11) tq => turb % uw
    if(i .eq. 12) tq => turb % vw
    tq % solver    = sstring
    tq % prec      = pstring
    tq % prec_opts = c_opts
    tq % tol       = tol
  end do

  end subroutine
