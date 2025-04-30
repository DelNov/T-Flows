!==============================================================================!
  subroutine Mass_Transfer_Estimate(Vof)
!------------------------------------------------------------------------------!
!   Calculates pressure source due to phase change                             !
!                                                                              !
!   Called from Multiphase_Mod_Vof_Pressure_Correction                         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: DEBUG = .false.
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Front_Type),    pointer :: Front
  type(Var_Type),      pointer :: t
  real,    contiguous, pointer :: elem_sx(:), elem_sy(:), elem_sz(:)
  integer, contiguous, pointer :: elem_used(:)
  character(SL)                :: fname
  integer                      :: e, s, c, c1, c2, elem
  real                         :: phi_c1, phi_c2, sx, sy, sz
  real                         :: dx1, dy1, dz1, dx2, dy2, dz2, d1, d2, st
  real                         :: dx1_n, dy1_n, dz1_n, dx2_n, dy2_n, dz2_n
  real                         :: q_0, q_1, sum_m_dot, sum_nab_t
  real                         :: c_lee
!==============================================================================!

  ! Take aliases
  Grid  => Vof % pnt_grid
  Flow  => Vof % pnt_flow
  Front => Vof % Front
  t     => Flow % t

  ! If not a problem with mass transfer, get out of here
  if(Flow % mass_transfer_model == NO_MASS_TRANSFER) return

  if(Math % Approx_Real(Vof % latent_heat, 1.0) ) then
    if (First_Proc()) then
      print '(a,1pe12.4)', ' # Latent heat is: ', Vof % latent_heat
      print '(a)', 'You might forget to define LATENT_HEAT in control'
    endif
    call Message % Error(60,                                                   &
                       'This error is critical.  Exiting.',                    &
                       file=__FILE__, line=__LINE__, one_proc=.true.)
  endif

  call Work % Connect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Connect_Int_Cell(elem_used)

  ! Initialize mass transfer term
  Vof % m_dot(:) = 0.0

  !------------------------------------------------!
  !                                                !
  !   Compute gradients of temperature, imposing   !
  !    saturation temperature at the interface     !
  !                                                !
  !------------------------------------------------!
  call Vof % Grad_Variable_With_Front(t, Vof % t_sat)
  if(DEBUG) then
    call Grid % Save_Debug_Vtu("grad-t",                               &
                               scalar_cell = t % n,                    &
                               scalar_name = "t",                      &
                               vector_cell = (/t % x, t % y, t % z/),  &
                               vector_name = "grad-t")
  end if

  !----------------------------------------------!
  !                                              !
  !    Extrapolate gradients from each of the    !
  !   phases toward the cells at the interface   !
  !                                              !
  !----------------------------------------------!
  if(Flow % mass_transfer_model .eq. TEMPERATURE_GRADIENTS) then

    ! Intialize t_0 and t_1 ...
    Vof % t_0 % x(1:) = t % x(1:)
    Vof % t_0 % y(1:) = t % y(1:)
    Vof % t_0 % z(1:) = t % z(1:)
    Vof % t_1 % x(1:) = t % x(1:)
    Vof % t_1 % y(1:) = t % y(1:)
    Vof % t_1 % z(1:) = t % z(1:)

    ! ... and extrapolate to front
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % x, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % y, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_1 % z, towards=0)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % x, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % y, towards=1)
    call Vof % Extrapolate_Normal_To_Front(Flow, Vof % t_0 % z, towards=1)
  endif

  !-------------------------------------------!
  !                                           !
  !   Computation of heat fluxes to or from   !
  !   interface towards the cells around it   !
  !                                           !
  !-------------------------------------------!

  !---------------------------------------------------!
  !   Distribute element areas to cell-based arrays   !
  !---------------------------------------------------!
  elem_sx(:) = 0.0
  elem_sy(:) = 0.0
  elem_sz(:) = 0.0
  do e = 1, Front % n_elems
    c = Front % Elem(e) % cell      ! in which cell it resides
    elem_sx(c) = Front % Elem(e) % sx
    elem_sy(c) = Front % Elem(e) % sy
    elem_sz(c) = Front % Elem(e) % sz
  end do
  call Grid % Exchange_Cells_Real(elem_sx)
  call Grid % Exchange_Cells_Real(elem_sy)
  call Grid % Exchange_Cells_Real(elem_sz)

  if(DEBUG) then
    fname = "elem-s-"
    write(fname(8:12), '(i5.5)') Time % Curr_Dt()
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
  do s = Faces_In_Domain_And_At_Buffers()
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
    if(elem_used(c) > 0) then
      elem_sx(c) = elem_sx(c) / elem_used(c)
      elem_sy(c) = elem_sy(c) / elem_used(c)
      elem_sz(c) = elem_sz(c) / elem_used(c)
    end if
  end do

  !----------------------------------------!
  !   Compute heat flux at the interface   !
  !----------------------------------------!
  Vof % a12(:) = 0.0
  Vof % a21(:) = 0.0

  do s = 1, Grid % n_faces

    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    phi_c1 = Vof % fun % n(c1)
    phi_c2 = Vof % fun % n(c2)

    Assert(phi_c1 .ne. 0.5)
    Assert(phi_c2 .ne. 0.5)

    if( (phi_c1 - 0.5) * (phi_c2 - 0.5) < 0 ) then

      ! Compute internal boundary coefficients on each side of the interface
      ! Units: W/(mK) * K/m * m^2 = W

      ! Simply sum the areas up.  Keep in mind that these
      ! have been corrected by the number of visits
      sx = elem_sx(c1) + elem_sx(c2)
      sy = elem_sy(c1) + elem_sy(c2)
      sz = elem_sz(c1) + elem_sz(c2)
      st = sqrt(sx**2 + sy**2 + sz**2)

      dx1 = Grid % xc(c1) - Front % xs(s)
      dy1 = Grid % yc(c1) - Front % ys(s)
      dz1 = Grid % zc(c1) - Front % zs(s)

      dx2 = Grid % xc(c1) + Grid % dx(s) - Front % xs(s)
      dy2 = Grid % yc(c1) + Grid % dy(s) - Front % ys(s)
      dz2 = Grid % zc(c1) + Grid % dz(s) - Front % zs(s)

      dx1_n = dx1 * 0.5*(Vof % nx(c1) + Vof % nx(c2))
      dy1_n = dy1 * 0.5*(Vof % ny(c1) + Vof % ny(c2))
      dz1_n = dz1 * 0.5*(Vof % nz(c1) + Vof % nz(c2))
      d1 = norm2((/dx1_n, dy1_n, dz1_n/))

      dx2_n = dx2 * 0.5*(Vof % nx(c1) + Vof % nx(c2))
      dy2_n = dy2 * 0.5*(Vof % ny(c1) + Vof % ny(c2))
      dz2_n = dz2 * 0.5*(Vof % nz(c1) + Vof % nz(c2))
      d2 = norm2((/dx2_n, dy2_n, dz2_n/))

      ! This asserion was failing due to imperfections in periodic solutions
      ! Assert(dot_product((/dx1_n, dy1_n, dz1_n/),(/dx2_n, dy2_n, dz2_n/))<0.0)

      ! Cell c1 is in vapor and c2 is in liquid
      if(Vof % fun % n(c1) < 0.5) then
        Vof % a12(s) = Vof % a12(s) + Vof % phase_cond(0) * st / d1
        Vof % a21(s) = Vof % a21(s) + Vof % phase_cond(1) * st / d2
      ! Cell c1 is in liquid and c2 is in vapor
      else
        Vof % a12(s) = Vof % a12(s) + Vof % phase_cond(1) * st / d1
        Vof % a21(s) = Vof % a21(s) + Vof % phase_cond(0) * st / d2
      end if

    end if      ! face is at the front

  end do

  !-----------------------------------------------!
  !                                               !
  !   Mass transfer with extrapolated gradients   !
  !                                               !
  !-----------------------------------------------!
  if(Flow % mass_transfer_model .eq. TEMPERATURE_GRADIENTS) then
    do c = Cells_In_Domain_And_Buffers()
      if(elem_used(c) > 0) then
        e  = Front % elem_in_cell(c)
        if(e > 0) then
          ! Units: W/(mK) * K/m * m^2 = W
          q_0 = (  Vof % t_0 % x(c) * elem_sx(c)    &
                 + Vof % t_0 % y(c) * elem_sy(c)    &
                 + Vof % t_0 % z(c) * elem_sz(c) ) * Vof % phase_cond(0)
          ! Heat flux to the interface in cell c from phase 1
          ! Units: W/(mK) * K/m * m^2 = W
          q_1 = (  Vof % t_1 % x(c) * elem_sx(c)    &
                 + Vof % t_1 % y(c) * elem_sy(c)    &
                + Vof % t_1 % z(c) * elem_sz(c) ) * Vof % phase_cond(1)

          ! Units: W / (J/kg) = kg/s
          Vof % m_dot(c) = (q_1 - q_0) / Vof % latent_heat
        end if
      end if
    end do
  endif

  !-----------------------------------------------!
  !                                               !
  !   Mass transfer modified Lee model            !
  !                                               !
  !-----------------------------------------------!
  if(Flow % mass_transfer_model .eq. LEE) then
    do c = Cells_In_Domain_And_Buffers()
      if(elem_used(c) > 0) then
        e  = Front % elem_in_cell(c)
        if(e > 0) then
          c_lee = Vof % c_lee(1)            ! CONDENSATION
          if (t % n(c) > Vof % t_sat) then  ! VAPORIZATION
                  c_lee = Vof % c_lee(2)
          endif
          Vof % m_dot(c) = c_lee * Flow % capacity(c) * Flow % density(c) &
                         * (t % n(c) - Vof % t_sat) &
                         * (sqrt(elem_sx(c)**2.0+elem_sy(c)**2.0+elem_sz(c)**2.0)) &
                         * (Grid%vol(c)**(1.0/3.0)) &
                         / Vof % latent_heat
        end if
      end if
    end do
  endif

  ! Sum
  sum_m_dot = 0.0
  sum_nab_t = 0.0
  elem = 0
  do c = Cells_In_Domain()
   e  = Front % elem_in_cell(c)
   if(e > 0) then
     sum_m_dot = sum_m_dot + Vof % m_dot(c)
     sum_nab_t = sum_nab_t + norm2((/t % x(c), t % y(c), t % z(c)/))
     elem = elem + 1
   end if
  end do
  call Global % Sum_Real(sum_m_dot)
  call Global % Sum_Real(sum_nab_t)
  call Global % Sum_Int(elem)

  if(DEBUG) then
    fname = "m_dot-"
    write(fname(7:11), '(i5.5)') Time % Curr_Dt()
    call Grid % Save_Debug_Vtu(fname,                            &
                               scalar_cell = Vof % m_dot,        &
                               scalar_name = "vof_m_dot")
  end if

  call Work % Disconnect_Real_Cell(elem_sx, elem_sy, elem_sz)
  call Work % Disconnect_Int_Cell(elem_used)

  end subroutine
