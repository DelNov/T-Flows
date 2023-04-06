!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof, curr_dt, ini)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  integer, intent(in)     :: curr_dt
  integer, intent(in)     :: ini
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
  real,    parameter :: L     = 2.26e+6  ! Latent heat [J/kg]   
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Var_Type),      pointer :: t
  real,    contiguous, pointer :: elem_sx(:), elem_sy(:), elem_sz(:)
  integer, contiguous, pointer :: elem_used(:)
  character(SL)                :: fname
  integer                      :: e, s, c, c1, c2
  real                         :: cond_1, cond_2, phi_c1, phi_c2, sx, sy, sz
!==============================================================================!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow
  t    => Flow % t

  call Work % Connect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Connect_Int_Cell(elem_used)

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  !------------------------------------------------!
  !   Compute gradients of temperature, imposing   !
  !    saturation temperature at the interface     !
  !------------------------------------------------!
  call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  if(DEBUG) then
    call Grid % Save_Debug_Vtu("grad-t",                               &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
  end if

  !---------------------------------------------------!
  !   Distribute element areas to cell-based arrays   !
  !---------------------------------------------------!
  elem_sx(:) = 0.0
  elem_sy(:) = 0.0
  elem_sz(:) = 0.0
  do e = 1, Vof % Front % n_elems
    c = Vof % Front % Elem(e) % cell      ! in which cell it resides
    elem_sx(c) = Vof % Front % Elem(e) % sx
    elem_sy(c) = Vof % Front % Elem(e) % sy
    elem_sz(c) = Vof % Front % Elem(e) % sz
  end do
  call Grid % Exchange_Cells_Real(elem_sx)
  call Grid % Exchange_Cells_Real(elem_sy)
  call Grid % Exchange_Cells_Real(elem_sz)
  if(DEBUG) then
    fname = "elem-s-"
    write(fname(8:12), '(i5.5)') curr_dt
    call Grid % Save_Debug_Vtu(fname,                      &
                               vector_cell = (/elem_sx,    &
                                               elem_sy,    &
                                               elem_sz/),  &
                               vector_name = "elem_s")
  end if

  !----------------------------------------------------!
  !   Count how many times will each element be used   !
  !----------------------------------------------------!
  elem_used(:) = 0
  do s = Faces_In_Domain()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_c1 = Vof % fun % n(c1)
    phi_c2 = Vof % fun % n(c2)

    ! If face crosses the 0.5 value increase the
    ! number of times the element has been used
    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then
      elem_used(c1) = elem_used(c1) + 1
      elem_used(c2) = elem_used(c2) + 1
    end if
  end do
  call Grid % Exchange_Cells_Int(elem_used)

  !------------------------------------------------------!
  !   Divide elements' surfaces by the number of times   !
  !     they are used to avoid excessive integration     !
  !------------------------------------------------------!
  do c = Cells_In_Domain_And_Buffers()
    elem_sx(c) = elem_sx(c) / elem_used(c)
    elem_sy(c) = elem_sy(c) / elem_used(c)
    elem_sz(c) = elem_sz(c) / elem_used(c)
  end do

  !----------------------------------------!
  !   Compute heat flux at the interface   !
  !----------------------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_c1 = Vof % fun % n(c1)
    phi_c2 = Vof % fun % n(c2)

    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then

      Vof % q_int(1,s) = 0.0
      Vof % q_int(2,s) = 0.0

      ! Take conductivities from each side of the interface
      cond_1 = Vof % phase_cond(1)
      cond_2 = Vof % phase_cond(0)
      if(Vof % fun % n(c1) < 0.5) cond_1 = Vof % phase_cond(0)
      if(Vof % fun % n(c2) > 0.5) cond_2 = Vof % phase_cond(1)

      ! Compute heat fluxes from each side of the interface
      ! Keep in mind that in the Stefan's case, dot product of
      ! element's surface and face's surface is positive.  If
      ! VOF's definition changes, I assume this would have to
      ! be adjusted accordingly.
      ! Units: W/(mK) * K/m * m^2 = W

      ! Simply sum the areas up.  Keep in mind that these
      ! have been corrected by the number of visits
      sx = elem_sx(c1) + elem_sx(c2)
      sy = elem_sy(c1) + elem_sy(c2)
      sz = elem_sz(c1) + elem_sz(c2)

      Vof % q_int(1,s) = Vof % q_int(1,s)  &
          + cond_1 * (t % x(c1) * sx + t % y(c1) * sy + t % z(c1) * sz)

      Vof % q_int(2,s) = Vof % q_int(2,s)  &
          + cond_2 * (t % x(c2) * sx + t % y(c2) * sy + t % z(c2) * sz)

    end if      ! face is at the front

  end do

  !-----------------------------------------------!
  !   Mass transfer with gradient extrapolation   !
  !-----------------------------------------------!

  ! Intialize t_0 and t_1 ...
  Vof % t_0 % x(1:) = t % x(1:)
  Vof % t_0 % y(1:) = t % y(1:)
  Vof % t_0 % z(1:) = t % z(1:)
  Vof % t_1 % x(1:) = t % x(1:)
  Vof % t_1 % y(1:) = t % y(1:)
  Vof % t_1 % z(1:) = t % z(1:)

  ! ... then extrapolate in the direction of the normal
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % x, towards=1)
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % y, towards=1)
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % z, towards=1)
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % x, towards=0)
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % y, towards=0)
  call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % z, towards=0)

  if(DEBUG) then
    fname = "grad-t-0-"
    write(fname(10:14), '(i5.5)') curr_dt
    call Grid % Save_Debug_Vtu(fname,                            &
                               scalar_cell = Vof % fun % n,      &
                               scalar_name = "vof_fun",          &
                               vector_cell = (/Vof % t_0 % x,    &
                                               Vof % t_0 % y,    &
                                               Vof % t_0 % z/),  &
                               vector_name = "grad-t-0")
    fname = "grad-t-1-"
    write(fname(10:14), '(i5.5)') curr_dt
    call Grid % Save_Debug_Vtu(fname,                            &
                               scalar_cell = Vof % fun % n,      &
                               scalar_name = "vof_fun",          &
                               vector_cell = (/Vof % t_1 % x,    &
                                               Vof % t_1 % y,    &
                                               Vof % t_1 % z/),  &
                               vector_name = "grad-t-1")
  end if

  Vof % q_0(:) = 0
  Vof % q_1(:) = 0
  do c = Cells_In_Domain_And_Buffers()

    if(elem_used(c) > 0) then
      ! Heat flux to the interface in cell c from phase 0
      ! Units: W/(mK) * K/m * m^2 = W
      Vof % q_0(c) = (  Vof % t_0 % x(c) * elem_sx(c)    &
                      + Vof % t_0 % y(c) * elem_sy(c)    &
                      + Vof % t_0 % z(c) * elem_sz(c) )  &
                   * Vof % phase_cond(0)

      ! Heat flux to the interface in cell c from phase 1
      ! Units: W/(mK) * K/m * m^2 = W
      Vof % q_1(c) = (  Vof % t_1 % x(c) * elem_sx(c)    &
                      + Vof % t_1 % y(c) * elem_sy(c)    &
                      + Vof % t_1 % z(c) * elem_sz(c) )  &
                   * Vof % phase_cond(1)

      ! Units: W / (J/kg) = kg/s
      Vof % m_dot(c) = (Vof % q_1(c) - Vof % q_0(c)) / L
    end if

  end do

  if(DEBUG) then
    fname = "m_dot-"
    write(fname(7:11), '(i5.5)') curr_dt
    call Grid % Save_Debug_Vtu(fname,                            &
                               scalar_cell = Vof % m_dot,        &
                               scalar_name = "vof_m_dot")
  end if

  call Work % Disconnect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Disconnect_Int_Cell(elem_used)

  end subroutine
