!==============================================================================!
  subroutine Multiphase_Mod_Compute_Vof(mult, sol, dt, n)
!------------------------------------------------------------------------------!
!   Solves Volume Fraction equation using UPWIND ADVECTION and CICSAM          !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: beta_f => r_face_01,  &
                      beta_c => r_face_02,  &
                      c_d    => r_cell_30
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
  integer                    :: s, c, c1, c2, fu
  integer                    :: donor, accept, corr_num
  integer                    :: i_sub, n_sub, wrong_vf, n_wrong_vf0, n_wrong_vf1
  integer,           pointer :: n_sub_param, corr_num_max
  character(len=80)          :: solver
  real                       :: fs
  real                       :: courant_max, epsloc
  real                       :: u_f, v_f, w_f
  real,              pointer :: courant_max_param
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

    if (vof % adv_scheme .eq. CICSAM) then
      ! Compute courant Number close to the interface:
      call Vof_Max_Courant_Number(mult, dt, c_d, 1, courant_max)

      n_sub = min(max(ceiling(courant_max / courant_max_param),1),n_sub_param)

      ! Warning if Courant Number is exceeded
      if (n_sub > 1) then
        call File_Mod_Append_File_For_Writing('alert-dt-vof.dat', fu)
        write(fu,*) 'Courant Number was exceded at iteration: ', n
        write(fu,*) 'Co_max = ', courant_max
        write(fu,*) 'Try reducing time step'
        close(fu)
      end if
    else
      n_sub = 1
    endif
  else
    ! Old volume fraction:
    vof % o(:) = vof % n(:)
  end if

  ! Phase change
  if (mult % phase_change) then
    call Multiphase_Mod_Vof_Mass_Transfer(mult, sol)
  end if

  if (vof % adv_scheme .eq. UPWIND) then
    !-------------------------!
    !   Matrix Coefficients   !
    !-------------------------!

    call Multiphase_Mod_Vof_Coefficients(mult, a, b, dt, beta_f)

    ! Solve System
    call Multiphase_Mod_Vof_Solve_System(mult, sol, b)

    call Grid_Mod_Exchange_Cells_Real(grid, vof % n)

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

      !---------------------------!
      !   Predict Beta at faces   !
      !---------------------------!

      ! Impose zero gradient at boundaries
      !call Multiphase_Mod_Vof_Boundary_Extrapolation(grid, mult, vof % n)
      do s = 1, grid % n_bnd_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
          vof % n(c2) = vof % n(c1)
        end if
      end do

      ! Old volume fraction:
      vof % o(:) = vof % n(:)

      ! Compute Gradient:
      call Field_Mod_Grad_Variable(flow, vof)

      ! Skewness correction
      call Multiphase_Mod_Vof_Skewness_Correction(mult, grid, beta_c)

      call Multiphase_Mod_Vof_Predict_Beta(mult, grid, beta_f, beta_c, c_d)

      loop_corr:  do corr_num = 1, corr_num_max
        !-------------------------!
        !   Matrix Coefficients   !
        !-------------------------!

        call Multiphase_Mod_Vof_Coefficients(mult, a, b,         &
                                             dt / real(n_sub),   &
                                             beta_f)

        ! Solve System
        call Multiphase_Mod_Vof_Solve_System(mult, sol, b)

        do s = 1, grid % n_bnd_faces
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
            vof % n(c2) = vof % n(c1)
          end if
        end do

        n_wrong_vf0 = 0
        n_wrong_vf1 = 0
        wrong_vf = 0

        ! Determine if 0 <= vof <= 1.0
        do c = 1, grid % n_cells
          if (vof % n(c) < -vof % tol) then
            n_wrong_vf0 = n_wrong_vf0 + 1
          end if

          if (vof % n(c) - 1.0 > vof % tol) then
            n_wrong_vf1 = n_wrong_vf1 + 1
          end if
        end do

        call Grid_Mod_Exchange_Cells_Real(grid, vof % n)

        !---------------------------!
        !   Correct Beta at faces   !
        !---------------------------!

        if (n_wrong_vf0 > 0 .or. n_wrong_vf1 > 0) then
          wrong_vf = 1
        end if

        call Comm_Mod_Global_Sum_Int(wrong_vf)

        if (wrong_vf == 0) then
          exit loop_corr
        else
          call Multiphase_Mod_Vof_Correct_Beta(mult, grid, beta_f, c_d)
        end if

      end do loop_corr

      !------------------------!
      !   Correct Boundaries   !
      !------------------------!
      do s = 1, grid % n_bnd_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
          if (vof % n(c2) < epsloc) then
            vof % n(c2) = 0.0
          end if
          if (vof % n(c2) - 1.0 >= epsloc) then
            vof % n(c2) = 1.0
          end if
        end if

      end do

      !--------------------------------------!
      !   Correct Interior Volume Fraction   !
      !--------------------------------------!
      do c = 1, grid % n_cells
        if (vof % n(c) < epsloc) then
          vof % n(c) = 0.0
        end if
        if (vof % n(c) - 1.0 >= epsloc) then
          vof % n(c) = 1.0
        end if
      end do

      call Grid_Mod_Exchange_Cells_Real(grid, vof % n)
    end do


  end if

  !------------------------------!
  !   Volume fraction at faces   !
  !------------------------------!

  if (vof % adv_scheme .eq. UPWIND) then

    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. INFLOW) then
        vof_f(s) = vof % n(c1)
        vof % n(c2) = vof % n(c1)
      end if
    end do

    ! Inside domain
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if (m_flux % n(s) >= 0.0) then
        vof_f(s) = vof % n(c1)
      else
        vof_f(s) = vof % n(c2)
      end if
    end do

  else if (vof % adv_scheme .eq. CICSAM .or. vof % adv_scheme .eq. STACS) then

    !call Multiphase_Mod_Vof_Boundary_Extrapolation(grid, mult, vof % n)
    ! At boundaries
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
        vof_f(s) = vof % n(c1)
        vof % n(c2) = vof % n(c1)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then
        vof_f(s) = vof % n(c1)
        vof % n(c2) = vof % n(c1)
      else if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        vof_f(s) = vof % n(c2)
      else
        vof_f(s) = vof % n(c1)
        vof % n(c2) = vof % n(c1)
      end if
    end do

    ! Inside domain
    do s = grid % n_bnd_faces + 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(m_flux % n(s) >= 0.0) then
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
    end do

  end if

  call Field_Mod_Grad_Variable(flow, vof)

  ! If distance function is calculated
  if (mult % d_func) then
    call Multiphase_Mod_Compute_Distance_Function(mult, sol, flow % dt, n)
  end if

  !----------------------------------------!
  !   Surface Tension Force Contribution   !
  !----------------------------------------!

  if (mult % surface_tension > TINY ) then
    if (mult % nodal_curvature) then
      call Multiphase_Mod_Vof_Surface_Tension_Contribution_Nodal(mult)
    else
      call Multiphase_Mod_Vof_Surface_Tension_Contribution_Csf(mult)
    end if
  end if

  !-----------------------!
  !   Update properties   !
  !-----------------------!
  call Multiphase_Mod_Update_Physical_Properties(mult, .false.)

  call Cpu_Timer_Mod_Stop('Compute_Multiphase (without solvers)')

  end subroutine
