!==============================================================================!
  subroutine Petsc_Solvers(Rc, Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   Reads options for PETSc solver from control file.                          !
!   I don't put it together with Read_Control_Numerical because I am afraid    !
!   it could grow a lot as we learn new PETSc options and want more from it    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc
  type(Field_Type),  target             :: Flow
  type(Turb_Type),   target             :: Turb
  type(Vof_Type),    target             :: Vof
  type(Solver_Type), target             :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, p, fun
  type(Var_Type),  pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: scalar(:)
  type(Var_Type),  pointer :: tq, phi
  logical                  :: found
  character(SL)            :: sstring, pstring
  character(SL)            :: opts(MSI)
  integer                  :: i, n_opts, sc
  real                     :: tol
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! If it wasn't compiled with PETSc, don't confuse a user with this
  if(.not. PETSC_ACTIVE) return

  ! Also, don't bother to read if user wants to use native solvers
  if(Sol % solvers .ne. PETSC) return

  ! Take aliases
  Grid   => Flow % pnt_grid
  t      => Flow % t
  p      => Flow % p
  scalar => Flow % scalar
  vis    => Turb % vis
  t2     => Turb % t2
  fun    => Vof % fun
  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)


  !----------------------------!
  !   For momentum equations   !
  !----------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_MOMENTUM',  &
                                        found,                        &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)

    u % solver = sstring
    v % solver = sstring
    w % solver = sstring
    u % prec = pstring
    v % prec = pstring
    w % prec = pstring
    u % prec_opts(1:MSI) = ''
    v % prec_opts(1:MSI) = ''
    w % prec_opts(1:MSI) = ''
    u % prec_opts(1:n_opts) = opts(1:n_opts)
    v % prec_opts(1:n_opts) = opts(1:n_opts)
    w % prec_opts(1:n_opts) = opts(1:n_opts)
    u % tol = tol
    v % tol = tol
    w % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for momentum are not specified.'  //  &
                ' Using the default values'
    end if
    u % prec = 'asm'
    v % prec = 'asm'
    w % prec = 'asm'
    u % prec_opts(1:MSI) = ''
    v % prec_opts(1:MSI) = ''
    w % prec_opts(1:MSI) = ''
  end if

  !---------------------------!
  !   For pressure equation   !
  !---------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_PRESSURE',  &
                                        found,                        &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'cg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm', pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol, .true.)

    Flow % pp % solver = sstring
    Flow % pp % prec = pstring
    Flow % pp % prec_opts(1:MSI)    = ''
    Flow % pp % prec_opts(1:n_opts) = opts(1:n_opts)
    Flow % pp % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for pressure are not specified.'  //  &
               ' Using the default values'
    end if
    Flow % pp % prec = 'asm'
    Flow % pp % prec_opts(1:MSI) = ''
  end if

  !-----------------------------------!
  !   For wall distance computation   !
  !-----------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_WALL_DISTANCE',  &
                                        found,                             &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol, .true.)

    Flow % wall_dist % solver = sstring
    Flow % wall_dist % prec = pstring
    Flow % wall_dist % prec_opts(1:MSI)    = ''
    Flow % wall_dist % prec_opts(1:n_opts) = opts(1:n_opts)
    Flow % wall_dist % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for potential are not specified.'  //  &
               ' Using the default values'
    end if
    Flow % wall_dist % prec = 'asm'
    Flow % wall_dist % prec_opts(1:MSI) = ''
  end if

  !----------------------------!
  !   For potential equation   !
  !----------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_POTENTIAL',  &
                                        found,                         &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol, .true.)

    Flow % pot % solver = sstring
    Flow % pot % prec = pstring
    Flow % pot % prec_opts(1:MSI)    = ''
    Flow % pot % prec_opts(1:n_opts) = opts(1:n_opts)
    Flow % pot % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for potential are not specified.'  //  &
               ' Using the default values'
    end if
    Flow % pot % prec = 'asm'
    Flow % pot % prec_opts(1:MSI) = ''
  end if

  !----------------------!
  !   For VOF function   !
  !----------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_VOF',  &
                                        found,                   &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol, .true.)

    Vof % fun % solver = sstring
    Vof % fun % prec = pstring
    Vof % fun % prec_opts(1:MSI)    = ''
    Vof % fun % prec_opts(1:n_opts) = opts(1:n_opts)
    Vof % fun % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for VOF are not specified.'  //  &
               ' Using the default values'
    end if
    Vof % fun % prec = 'asm'
    Vof % fun % prec_opts(1:MSI) = ''
  end if

  !-------------------------!
  !   For energy equation   !
  !-------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_ENERGY',  &
                                        found,                      &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)

    t % solver = sstring
    t % prec = pstring
    t % prec_opts(1:MSI)    = ''
    t % prec_opts(1:n_opts) = opts(1:n_opts)
    t % tol = tol
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for energy are not specified.'  //  &
               ' Using the default values'
    end if
    t % prec = 'asm'
    t % prec_opts(1:MSI) = ''
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_SCALARS',  &
                                        found,                       &
                                        .false.)

  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)

    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      phi % solver = sstring
      phi % prec = pstring
      phi % prec_opts(1:MSI)    = ''
      phi % prec_opts(1:n_opts) = opts(1:n_opts)
      phi % tol = tol
    end do
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for scalars are not specified.'  //  &
               ' Using the default values'
    end if
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      phi % prec = 'asm'
      phi % prec_opts(1:MSI) = ''
    end do
  end if

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_TURBULENCE',  &
                                        found,                          &
                                        .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER', 'bicg', sstring, .true.)
    call Control % Read_Char_Item_On('PREC',   'asm',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', opts, n_opts, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol, .true.)

    do i = 1, 12
      if(i .eq.  1) tq => Turb % kin
      if(i .eq.  2) tq => Turb % eps
      if(i .eq.  3) tq => Turb % zeta
      if(i .eq.  4) tq => Turb % f22
      if(i .eq.  5) tq => Turb % vis
      if(i .eq.  6) tq => Turb % t2
      if(i .eq.  7) tq => Turb % uu
      if(i .eq.  8) tq => Turb % vv
      if(i .eq.  9) tq => Turb % ww
      if(i .eq. 10) tq => Turb % uv
      if(i .eq. 11) tq => Turb % uw
      if(i .eq. 12) tq => Turb % vw
      tq % solver              = sstring
      tq % prec                = pstring
      tq % prec_opts(1:MSI)    = ''
      tq % prec_opts(1:n_opts) = opts(1:n_opts)
      tq % tol                 = tol
    end do
  else
    if(First_Proc()) then
      print *, '# NOTE! PETSc options for turbulence are not specified.'  //  &
               ' Using the default values'
    end if
    do i = 1, 12
      if(i .eq.  1) tq => Turb % kin
      if(i .eq.  2) tq => Turb % eps
      if(i .eq.  3) tq => Turb % zeta
      if(i .eq.  4) tq => Turb % f22
      if(i .eq.  5) tq => Turb % vis
      if(i .eq.  6) tq => Turb % t2
      if(i .eq.  7) tq => Turb % uu
      if(i .eq.  8) tq => Turb % vv
      if(i .eq.  9) tq => Turb % ww
      if(i .eq. 10) tq => Turb % uv
      if(i .eq. 11) tq => Turb % uw
      if(i .eq. 12) tq => Turb % vw
      tq % prec = 'asm'
      tq % prec_opts(1:MSI) = ''
    end do
  end if

  end subroutine
