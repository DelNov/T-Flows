!==============================================================================!
  subroutine Multiphase_Mod_Compute_Vof(mult, sol, dt, n)
!------------------------------------------------------------------------------!
!   Solves Volume Fraction equation using UPWIND ADVECTION and CICSAM          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: beta_f => r_face_01, c_d => r_cell_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  real                          :: dt
  integer                       :: n    ! current temporal iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: vof
  type(Face_Type),   pointer :: m_flux
  real, contiguous,  pointer :: vof_f(:)
  type(Matrix_Type), pointer :: a
  real, contiguous,  pointer :: b(:)
  integer                    :: s, c, c1, c2
  integer                    :: donor, accept, corr_num
  integer                    :: i_sub, n_sub
  real                       :: fs
  real                       :: courant_max, epsloc
  real,              pointer :: courant_max_param
  integer,           pointer :: n_sub_param, corr_num_max
  character(len=80)          :: solver
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Multiphase (without solvers)')

  ! Take aliases
  flow   => mult % pnt_flow
  grid   => flow % pnt_grid
  m_flux => flow % m_flux
  vof    => mult % vof
  vof_f  => mult % vof_f
  courant_max_param => mult % courant_max_param
  n_sub_param       => mult % n_sub_param
  corr_num_max      => mult % corr_num_max

  a => sol % a
  b => sol % b % val

  epsloc = epsilon(epsloc)

  if (vof % adv_scheme .eq. CICSAM .or. &
      vof % adv_scheme .eq. STACS) then

    ! Compute courant Number close to the interface:
    call Vof_Max_Courant_Number(mult, dt, c_d, 1, courant_max)

    n_sub = min(max(ceiling(courant_max / courant_max_param),1),n_sub_param)

  else
    ! Old volume fraction:
    vof % o(:) = vof % n(:)
  end if

  if (vof % adv_scheme .eq. UPWIND) then
    !-------------------------!
    !   Matrix Coefficients   !
    !-------------------------!

    call Multiphase_Mod_Vof_Coefficients(mult, a, b, dt, beta_f)

    ! Solve System
    call Multiphase_Mod_Vof_Solve_System(mult, sol, b)

    call Grid_Mod_Exchange_Real(grid, vof % n)

    !-----------------------------!
    !   Correct Volume Fraction   !
    !-----------------------------!

    do c = 1, grid % n_cells
      vof % n(c) = max(min(vof % n(c),1.0),0.0)
    end do

  else if (vof % adv_scheme .eq. CICSAM .or. &
           vof % adv_scheme .eq. STACS) then

    do i_sub = 1, n_sub

      ! Courant number full domain:
      call Vof_Max_Courant_Number(mult, dt / real(n_sub),    &
                                  c_d, 0, courant_max)

      ! Old volume fraction:
      vof % o(:) = vof % n(:)

      !---------------------------!
      !   Predict Beta at faces   !
      !---------------------------!
      ! Compute Gradient:
      call Field_Mod_Grad_Variable(flow, vof)

      call Multiphase_Mod_Vof_Predict_Beta(vof,                        &
                                           m_flux % n,                 &
                                           vof % x, vof % y, vof % z,  &
                                           grid % dx,                  &
                                           grid % dy,                  &
                                           grid % dz,                  &
                                           beta_f,                     &
                                           c_d)

      loop_corr:  do corr_num = 1, corr_num_max
        !-------------------------!
        !   Matrix Coefficients   !
        !-------------------------!

        call Multiphase_Mod_Vof_Coefficients(mult, a, b,         &
                                             dt / real(n_sub),   &
                                             beta_f)

        ! Solve System
        call Multiphase_Mod_Vof_Solve_System(mult, sol, b)

        call Grid_Mod_Exchange_Real(grid, vof % n)

        !---------------------------!
        !   Correct Beta at faces   !
        !---------------------------!
        call Multiphase_Mod_Vof_Correct_Beta(mult,               &
                                             vof,                &
                                             m_flux % n,         &
                                             beta_f,             &
                                             dt / real(n_sub))

        !------------------------!
        !   Correct Boundaries   !
        !------------------------!
        do s = 1, grid % n_faces

          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)

          if (c2 < 0) then
            if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
              vof % n(c2) = max(min(vof % n(c1),1.0),0.0)
            end if
            if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL) then
              vof % n(c2) = max(min(vof % n(c1),1.0),0.0)
            end if
          end if

        end do

        !---------------------!
        !   Noise Reduction   !
        !---------------------!
        do c = 1, grid % n_cells
          if(vof % n(c) < epsloc) then
            vof % n(c) = 0.0
          end if

          if(1.0 - vof % n(c) < epsloc) then
            vof % n(c) = 1.0
          end if
        end do

      end do loop_corr

      !-----------------------------!
      !   Correct Volume Fraction   !
      !-----------------------------!
      do c = 1, grid % n_cells
        vof % n(c) = max(min(vof % n(c),1.0),0.0)
      end do

      call Grid_Mod_Exchange_Real(grid, vof % n)
    end do


  end if

  !------------------------------!
  !   Volume fraction at faces   !
  !------------------------------!

  if (vof % adv_scheme .eq. UPWIND) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      ! Face is inside the domain
      if(c2 > 0) then

        if (m_flux % n(s)>=0.0) then
          vof_f(s) = vof % n(c1)
        else
          vof_f(s) = vof % n(c2)
        end if

      ! Side is on the boundary
      else ! (c2 < 0)

        vof_f(s) = vof % n(c1)

      end if

    end do

  else if (vof % adv_scheme .eq. CICSAM .or. vof % adv_scheme .eq. STACS) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      fs = grid % f(s)

      if(c2 > 0) then
        if(m_flux % n(s) >= 0) then
          donor = c1
          accept = c2
        else
          donor = c2
          accept = c1
        end if

        vof_f(s) = 0.5 * ((1.0 - beta_f(s)) * ( vof % n(donor)      &
                                              + vof % o(donor) )    &
                               + beta_f(s)  * ( vof % n(accept)     &
                                              + vof % o(accept) ))

      else
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          vof_f(s) = vof % n(c1)
        else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          vof_f(s) = vof % n(c2)
        else
          vof_f(s) = vof % n(c1)
        end if

      end if
    end do

  end if

  !-----------------------!
  !   Update properties   !
  !-----------------------!
  call Multiphase_Mod_Update_Physical_Properties(mult)

  call Field_Mod_Grad_Variable(flow, vof)

  !----------------------------------------!
  !   Surface Tension Force Contribution   !
  !----------------------------------------!

  ! If distance function is calculated
  if (mult % d_func) then
    call Multiphase_Mod_Compute_Distance_Function(mult, sol, flow % dt, n)
  end if

  if (mult % surface_tension > TINY ) then
    call Multiphase_Mod_Vof_Surface_Tension_Contribution(mult)
  end if

  call Cpu_Timer_Mod_Stop('Compute_Multiphase (without solvers)')

  end subroutine
