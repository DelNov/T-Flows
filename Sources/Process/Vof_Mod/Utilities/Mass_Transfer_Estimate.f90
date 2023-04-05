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
  logical, parameter :: DEBUG = .true.
  logical, parameter :: EXTRAPOLATION = .true.
  real,    parameter :: L = 2.26e+6  ! Latent heat [J/kg]   
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer     :: Grid
  type(Field_Type), pointer     :: Flow
  type(Var_Type),   pointer     :: t
  integer                       :: e, s, c1, c2, i_ele
  real                          :: cond_1, cond_2, tmp
  integer                       :: c
  character(SL)                 :: fname
  integer                       :: n, cn, i_cel, i_nod, i_fac, j_cel
  real                          :: x, y, z, xe, ye, ze, xcn, ycn, zcn
  real                          :: dx, dy, dz, c_size, dg
!==============================================================================!

  !------------------------------------------------!
  !   MASS TRANSFER OLD PART WITH EXTRAPOLATION    !
  !------------------------------------------------!

  ! Take aliases
  Grid => Vof % pnt_grid
  Flow => Vof % pnt_flow
  t    => Flow % t

  ! If not a problem with mass transfer, get out of here
  if(.not. Flow % mass_transfer) return

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  !------------------------------------------------!
  !  Compute gradients of temperature, imposing   !
  !   saturation temperature at the interface     !
  !------------------------------------------------!
  call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  
  if(DEBUG) then
    call Grid % Save_Debug_Vtu("grad-t",                               &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
  end if

  !----------------------------------------!
  !  Compute heat flux at the interface   !
  !----------------------------------------!
  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    if(any(Vof % Front % elems_at_face(1:2,s) .ne. 0)) then

      Vof % q_int(1,s) = 0.0
      Vof % q_int(2,s) = 0.0

      do i_ele = 1, 2
        e = Vof % Front % elems_at_face(i_ele, s)
        if(e .ne. 0) then

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
          if(  Vof % Front % elem(e) % sx * Grid % sx(s)  &
             + Vof % Front % elem(e) % sy * Grid % sy(s)  &
             + Vof % Front % elem(e) % sz * Grid % sz(s) > 0.0 ) then

               Vof % q_int(1,s) = Vof % q_int(1,s)                        &
                    + cond_1 * ( t % x(c1) * Vof % Front % elem(e) % sx   &
                               + t % y(c1) * Vof % Front % elem(e) % sy   &
                               + t % z(c1) * Vof % Front % elem(e) % sz)

               Vof % q_int(2,s) = Vof % q_int(2,s)                        &
                    + cond_2 * ( t % x(c2) * Vof % Front % elem(e) % sx   &
                               + t % y(c2) * Vof % Front % elem(e) % sy   &
                               + t % z(c2) * Vof % Front % elem(e) % sz)
          else

               Vof % q_int(1,s) = Vof % q_int(1,s)                        &
                    - cond_1 * ( t % x(c1) * Vof % Front % elem(e) % sx   &
                               + t % y(c1) * Vof % Front % elem(e) % sy   &
                               + t % z(c1) * Vof % Front % elem(e) % sz)

               Vof % q_int(2,s) = Vof % q_int(2,s)                        &
                    - cond_2 * ( t % x(c2) * Vof % Front % elem(e) % sx   &
                               + t % y(c2) * Vof % Front % elem(e) % sy   &
                               + t % z(c2) * Vof % Front % elem(e) % sz)
          end if

        end if  ! e .ne. 0
      end do    ! i_ele

    end if      ! face is at the front

  end do

  !-------------------------------------------------!
  !   MASS TRANSFER WITH GRADIENT EXTRAPOLATION     !
  !-------------------------------------------------!
  if(EXTRAPOLATION) then
  ! Gradients are fresh here
  ! Intialize t_0 and t_1
    Vof % t_0 % x(1:) = Flow % t % x(1:)
    Vof % t_0 % y(1:) = Flow % t % y(1:)
    Vof % t_0 % z(1:) = Flow % t % z(1:)
    Vof % t_1 % x(1:) = Flow % t % x(1:)
    Vof % t_1 % y(1:) = Flow % t % y(1:)
    Vof % t_1 % z(1:) = Flow % t % z(1:)
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
    do c = Cells_In_Domain()  ! Cells_In_Domain should suffice
      e = Vof % Front % elem_in_cell(c)
      if(e > 0) then
        ! Heat flux to the interface in cell c from phase 0
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_0(c) = (  Vof % t_0 % x(c) * Vof % Front % Elem(e) % sx   &
                        + Vof % t_0 % y(c) * Vof % Front % Elem(e) % sy   &
                        + Vof % t_0 % z(c) * Vof % Front % Elem(e) % sz)  &
                     * Vof % phase_cond(0)

        ! Heat flux to the interface in cell c from phase 1
        ! Units: W/(mK) * K/m * m^2 = W
        Vof % q_1(c) = (  Vof % t_1 % x(c) * Vof % Front % Elem(e) % sx   &
                        + Vof % t_1 % y(c) * Vof % Front % Elem(e) % sy   &
                        + Vof % t_1 % z(c) * Vof % Front % Elem(e) % sz)  &
                     * Vof % phase_cond(1)
        ! Units: W / (J/kg) = kg/s
        Vof % m_dot(c) = (Vof % q_1(c) - Vof % q_0(c)) / L
      end if
    end do
    call Grid % Exchange_Cells_Real(Vof % q_0)
    call Grid % Exchange_Cells_Real(Vof % q_1)
    if(DEBUG) then
      fname = "m_dot-"
      write(fname(7:11), '(i5.5)') curr_dt
      call Grid % Save_Debug_Vtu(fname,                            &
                                 scalar_cell = Vof % m_dot,        &
                                 scalar_name = "vof_m_dot")
    end if
  end if

  end subroutine
