!==============================================================================!
  subroutine Petsc_Solvers(Rc, Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to read and set up the parameters for various
!>  PETSc solvers from the control file.  It configures solvers for different
!>  aspects of the flow simulation such as momentum, pressure, wall distance,
!>  potential, heat transfer, multiphase flow, passive scalars, turbulence
!>  quantities, and also reads options to be passed to PETSc library.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type), intent(in) :: Rc    !! parent class
  type(Field_Type),  target             :: Flow  !! flow object
  type(Turb_Type),   target             :: Turb  !! turbulence object
  type(Vof_Type),    target             :: Vof   !! VOF object
  type(Solver_Type), target             :: Sol   !! solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, p, fun
  type(Var_Type),  pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: scalar(:)
  type(Var_Type),  pointer :: tq, phi
  logical                  :: found
  character(SL)            :: sstring, pstring
  character(SL)            :: pets(MAX_STRING_ITEMS), prec(MAX_STRING_ITEMS)
  integer                  :: i, n_pets, n_prec, sc
  real                     :: tol
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! If it wasn't compiled with PETSc, don't confuse a user with this
  if(.not. PETSC_ACTIVE) return

  ! Give some sign
  if(First_Proc())  &
    print '(a)', ' # Reading info about PETSc solvers'

  ! Initialize PETSc options to be all empty strings
  petsc_options(1:MAX_STRING_ITEMS) = ''

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

  !---------------------------!
  !   General PETSc options   !
  !---------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS', found, .false.)
  if(found) then
    call Control % Read_Strings_On('VALUES', pets, n_pets, .false.)
    petsc_options(1:n_pets) = pets(1:n_pets)
  end if

  !----------------------------!
  !   For momentum equations   !
  !----------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_MOMENTUM',  &
                                      found,                        &
                                      .false.)
  u % o_prec(1:MAX_STRING_ITEMS) = ''
  v % o_prec(1:MAX_STRING_ITEMS) = ''
  w % o_prec(1:MAX_STRING_ITEMS) = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'bicg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'asm',   pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol,     .true.)

    u % solver = sstring
    v % solver = sstring
    w % solver = sstring
    u % prec = pstring
    v % prec = pstring
    w % prec = pstring
    u % o_prec(1:n_prec) = prec(1:n_prec)
    v % o_prec(1:n_prec) = prec(1:n_prec)
    w % o_prec(1:n_prec) = prec(1:n_prec)
    u % tol = tol
    v % tol = tol
    w % tol = tol
  else
    u % solver = 'bicg'
    v % solver = 'bicg'
    w % solver = 'bicg'
    u % prec   = 'asm'
    v % prec   = 'asm'
    w % prec   = 'asm'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for momentum are not'  //  &
                   ' specified.  Using the default: '             //  &
                   trim(u % solver) // '/' // trim(u % prec)
    end if
  end if

  !---------------------------!
  !   For pressure equation   !
  !---------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_PRESSURE',  &
                                      found,                        &
                                      .false.)
  Flow % pp % o_prec(1:MAX_STRING_ITEMS)    = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'cg',    sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'gamg',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol,     .true.)

    Flow % pp % solver = sstring
    Flow % pp % prec = pstring
    Flow % pp % o_prec(1:n_prec) = prec(1:n_prec)
    Flow % pp % tol = tol
  else
    Flow % pp % solver = 'cg'
    Flow % pp % prec   = 'gamg'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for pressure are not'  //  &
                   ' specified.  Using the default: '             //  &
                   trim(Flow % pp % solver) // '/' // trim(Flow % pp % prec)
    end if
  end if

  !-----------------------------------!
  !   For wall distance computation   !
  !-----------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_WALL_DISTANCE',  &
                                        found,                             &
                                        .false.)
  Flow % wall_dist % o_prec(1:MAX_STRING_ITEMS)    = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'cg',    sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'gamg',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol,     .true.)

    Flow % wall_dist % solver = sstring
    Flow % wall_dist % prec = pstring
    Flow % wall_dist % o_prec(1:n_prec) = prec(1:n_prec)
    Flow % wall_dist % tol = tol
  else
    Flow % wall_dist % solver = 'cg'
    Flow % wall_dist % prec   = 'gamg'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for wall distance are not'  //  &
                   ' specified.  Using the default: '                  //  &
                   trim(Flow % wall_dist % solver) // '/'              //  &
                   trim(Flow % wall_dist % prec)
    end if
  end if

  !----------------------------!
  !   For potential equation   !
  !----------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_POTENTIAL',  &
                                      found,                         &
                                      .false.)
  Flow % pot % o_prec(1:MAX_STRING_ITEMS) = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'cg',    sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'gamg',  pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol,     .true.)

    Flow % pot % solver = sstring
    Flow % pot % prec = pstring
    Flow % pot % o_prec(1:n_prec) = prec(1:n_prec)
    Flow % pot % tol = tol
  else
    Flow % pot % solver = 'cg'
    Flow % pot % prec   = 'gamg'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for potential are not'  //  &
                   ' specified.  Using the default: '              //  &
                   trim(Flow % pot % solver) // '/' // trim(Flow % pot % prec)
    end if
  end if

  !----------------------!
  !   For VOF function   !
  !----------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_VOF',  &
                                      found,                   &
                                      .false.)
  Vof % fun % o_prec(1:MAX_STRING_ITEMS) = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'bicg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'asm',   pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol,     .true.)

    Vof % fun % solver = sstring
    Vof % fun % prec = pstring
    Vof % fun % o_prec(1:n_prec) = prec(1:n_prec)
    Vof % fun % tol = tol
  else
    Vof % fun % solver = 'bicg'
    Vof % fun % prec   = 'asm'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for VOF are not'  //  &
                   ' specified.  Using the default: '        //  &
                   trim(Vof % fun % solver) // '/' // trim(Vof % fun % prec)
    end if
  end if

  !-------------------------!
  !   For energy equation   !
  !-------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_ENERGY',  &
                                      found,                      &
                                      .false.)
  t % o_prec(1:MAX_STRING_ITEMS) = ''
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'bicg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'asm',   pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-5, tol,     .true.)

    t % solver = sstring
    t % prec = pstring
    t % o_prec(1:n_prec) = prec(1:n_prec)
    t % tol = tol
  else
    t % solver = 'bicg'
    t % prec   = 'asm'
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for energy are not'  //  &
                   ' specified.  Using the default: '           //  &
                   trim(t % solver) // '/' // trim(t % prec)
    end if
  end if

  !--------------------------------!
  !   Related to passive scalars   !
  !--------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_SCALARS',  &
                                      found,                       &
                                      .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'bicg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'asm',   pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol,     .true.)

    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      phi % solver = sstring
      phi % prec = pstring
      phi % o_prec(1:n_prec) = prec(1:n_prec)
      phi % o_prec(1:MAX_STRING_ITEMS) = ''
      phi % tol = tol
    end do
  else
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      phi % solver = 'bicg'
      phi % prec   = 'asm'
      phi % o_prec(1:MAX_STRING_ITEMS) = ''
    end do
    if(Flow % n_scalars > 0) then
      phi => Flow % scalar(1)  ! probably not needed, but doesn't harm
      if(First_Proc() .and. associated(phi)) then
        print '(a)', ' # NOTE! PETSc options for scalars are not'  //  &
                     ' specified.  Using the default: '            //  &
                     trim(phi % solver) // '/' // trim(phi % prec)
      end if
    end if
  end if

  !------------------------------!
  !   All turbuelnt quantities   !
  !------------------------------!
  call Control % Position_At_One_Key('PETSC_OPTIONS_FOR_TURBULENCE',  &
                                      found,                          &
                                      .false.)
  if(found) then
    call Control % Read_Char_Item_On('SOLVER',   'bicg',  sstring, .true.)
    call Control % Read_Char_Item_On('PREC',     'asm',   pstring, .true.)
    call Control % Read_Strings_On  ('PREC_OPTS', prec,   n_prec, .false.)
    call Control % Read_Real_Item_On('TOLERANCE', 1.0e-3, tol,     .true.)

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
      tq % solver           = sstring
      tq % prec             = pstring
      tq % o_prec(:)        = ''
      tq % o_prec(1:n_prec) = prec(1:n_prec)
      tq % tol              = tol
    end do
  else
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
      tq % solver    = 'bicg'
      tq % prec      = 'asm'
      tq % o_prec(:) = ''
    end do
    if(First_Proc()) then
      print '(a)', ' # NOTE! PETSc options for turbulence are not'  //  &
                   ' specified.  Using the default: '               //  &
                   trim(tq % solver) // '/' // trim(tq % prec)
    end if
  end if

  end subroutine
