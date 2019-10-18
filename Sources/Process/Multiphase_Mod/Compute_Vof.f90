!======================================================================!
  subroutine Multiphase_Mod_Compute_Vof(mult, sol, dt, ini)
!----------------------------------------------------------------------!
!   Forms and solves Volume Fraction equation using UPWIND ADVECTION   !
!----------------------------------------------------------------------!
!------------------------------[Modules]-------------------------------!
  use Const_Mod
  use Comm_Mod
  use Cpu_Timer_Mod, only: Cpu_Timer_Mod_Start, Cpu_Timer_Mod_Stop
  use Field_Mod,     only: Field_Type, density, viscosity, dens_face
  use Var_Mod,       only: Var_Type
  use Grid_Mod,      only: Grid_Type
  use Bulk_Mod,      only: Bulk_Type
  use Info_Mod,      only: Info_Mod_Iter_Fill_At
  use Numerics_Mod
  use Solver_Mod,    only: Solver_Type, Bicg, Cg, Cgs, Acm
  use Matrix_Mod,    only: Matrix_Type
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Work_Mod, only: beta_f => r_face_01
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  real                          :: dt
  integer                       :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type),  pointer :: flow
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: vof
  type(Var_type),    pointer :: vol_flux
  real,              pointer :: vof_f(:)
  real,              pointer :: vof_i(:), vof_j(:), vof_k(:)
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: s, c, c1, c2
  integer                    :: exec_iter
  integer                    :: donor, accept, corr_num, corr_num_max
  real                       :: fs, a0
  character(len=80)          :: solver
  real                       :: upwd1, upwd2, beta_dummy
  logical                    :: corr_cicsam
!==============================================================================!

  call Cpu_Timer_Mod_Start('Compute_Multiphase (without solvers)')

  ! Take aliases
  flow     => mult % pnt_flow
  grid     => flow % pnt_grid
  bulk     => flow % bulk
  vof      => mult % vof
  vof_f    => mult % vof_f
  vol_flux => flow % vol_flux
  vof_i    => vof % x
  vof_j    => vof % y
  vof_k    => vof % z
  flux     => flow % flux

  a => sol % a
  b => sol % b % val

  if (vof % adv_scheme .eq. CICSAM .or. &
      vof % adv_scheme .eq. STACS) then
    ! Stablish max and min:
    call Numerics_Mod_Min_Max(vof)

    corr_num_max = 5
  end if

  ! Old Volume fraction and flux:
  if(ini .eq. 1) then
    vof % o(:) = vof % n(:)
  end if


  if (vof % adv_scheme .eq. UPWIND) then
    !-------------------------!
    !   Matrix Coefficients   !
    !-------------------------!

    ! Initialize matrix and right hand side
    b       = 0.0
    a % val = 0.0

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      upwd1 = 0.5 * max(vol_flux % n(s), 0.0)
      upwd2 = 0.5 * max(-vol_flux % n(s), 0.0)

      a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1
      b(c1) = b(c1) - ( upwd1 * vof % o(c1) -  upwd2 * vof % o(c2) ) 

      if (c2 > 0) then
        a % val(a % pos(1,s)) =  - upwd2

        a % val(a % dia(c2)) = a % val(a % dia(c2)) + upwd2
        b(c2) = b(c2) - ( upwd2 * vof % o(c2) - upwd1 * vof % o(c1))
        a % val(a % pos(2,s)) =  - upwd1 
      else
        if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          b(c1) = b(c1) - vol_flux % n(s) * vof % n(c2)
        else if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
          a % val(a % dia(c1)) = a % val(a % dia(c1)) + vol_flux % n(s)
        end if
      end if

    end do

    !------------------------------------------------------------!
    !                Calculate Main coefficient and              !
    !   Calculate Source term (part 2 : Temporal contribution)   !
    !------------------------------------------------------------!

    ! Two time levels; Linear interpolation
    if(vof % td_scheme .eq. LINEAR) then
      do c = 1, grid % n_cells
        a0 = grid % vol(c) / dt
        a % val(a % dia(c)) = a % val(a % dia(c)) + a0
        b(c)  = b(c) + a0 * vof % o(c)
      end do
    end if

    ! Three time levels; parabolic interpolation
    if(vof % td_scheme .eq. PARABOLIC) then
      do c = 1, grid % n_cells
        a0 = grid % vol(c) / dt
        a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
        b(c)  = b(c) + 2.0 * a0 * vof % o(c) - 0.5 * a0 * vof % oo(c)
      end do
    end if

    ! Get solver
    call Control_Mod_Solver_For_Multiphase(solver)

    call Cpu_Timer_Mod_Start('Linear_Solver_For_Multiphase')
    if(solver .eq. 'ACM') then
      vof % tol   = PICO
      call Acm(sol,            &
               vof % n,        &
               b,              &
               vof % precond,  &
               vof % niter,    &     ! number of V cycles
               vof % tol,      &
               vof % res)
      stop
    else
      call Cg(sol,            &
              vof % n,        &
              b,              &
              vof % precond,  &
              vof % niter,    &      ! max number of iterations
              exec_iter,      &      ! executed number of iterations
              vof % tol,      &
              vof % res)
    end if

    call Cpu_Timer_Mod_Stop('Linear_Solver_For_Multiphase')

    call Info_Mod_Iter_Fill_At(1, 6, vof % name, exec_iter, vof % res)

    !-----------------------------!
    !   Correct Volume Fraction   !
    !-----------------------------!

    do c = 1, grid % n_cells
      if(vof % n(c)>1.0) then
        vof % n(c) = 1.0
      end if

      if(vof % n(c)<0.0) then
        vof % n(c) = 0.0
      end if
    end do

  else if (vof % adv_scheme .eq. CICSAM .or. &
           vof % adv_scheme .eq. STACS) then
    !---------------------------!
    !   Predict Beta at faces   !
    !---------------------------!
    ! Compute Gradient:
    call Grad_Mod_Variable(vof)

    call Multiphase_Mod_Vof_Predict_Beta(vof,                  &
                                         vol_flux % n,         &
                                         vof_i, vof_j, vof_k,  &
                                         grid % dx,            &
                                         grid % dy,            &
                                         grid % dz,            &
                                         beta_f,               &
                                         dt)
  loop_corr:  do corr_num = 1, corr_num_max
      !-------------------------!
      !   Matrix Coefficients   !
      !-------------------------!

      ! Initialize matrix and right hand side
      b       = 0.0
      a % val = 0.0

      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        if (vol_flux % n(s) > 0.0) then
          beta_dummy = beta_f(s)
        else
          beta_dummy = 1.0 - beta_f(s)
        end if
        upwd1 = 0.5 * (1.0 - beta_dummy) * vol_flux % n(s)
        upwd2 = 0.5 * beta_dummy * vol_flux % n(s)

        if (c2 > 0) then
          a % val(a % dia(c1)) = a % val(a % dia(c1)) + upwd1
          a % val(a % pos(1,s)) = upwd2
          b(c1) = b(c1) - upwd1 * vof % o(c1) - upwd2 * vof % o(c2)

          a % val(a % dia(c2)) = a % val(a % dia(c2)) - upwd2 
          a % val(a % pos(2,s)) =  - upwd1 
          b(c2) = b(c2) + upwd2 * vof % o(c2) + upwd1 * vof % o(c1)
        else
          if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
            b(c1) = b(c1) - vol_flux % n(s) * vof % n(c2)
          else if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
            a % val(a % dia(c1)) = a % val(a % dia(c1)) + vol_flux % n(s)
          end if
        end if

      end do

      !------------------------------------------------------------!
      !                Calculate Main coefficient and              !
      !   Calculate Source term (part 2 : Temporal contribution)   !
      !------------------------------------------------------------!

      ! Two time levels; Linear interpolation
      if(vof % td_scheme .eq. LINEAR) then
        do c = 1, grid % n_cells
          a0 = grid % vol(c) / dt
          a % val(a % dia(c)) = a % val(a % dia(c)) + a0
          b(c)  = b(c) + a0 * vof % o(c)
        end do
      end if

      ! Three time levels; parabolic interpolation
      if(vof % td_scheme .eq. PARABOLIC) then
        do c = 1, grid % n_cells
          a0 = grid % vol(c) / dt
          a % val(a % dia(c)) = a % val(a % dia(c)) + 1.5 * a0
          b(c)  = b(c) + 2.0 * a0 * vof % o(c) - 0.5 * a0 * vof % oo(c)
        end do
      end if

      ! Get solver
      call Control_Mod_Solver_For_Multiphase(solver)

      call Cpu_Timer_Mod_Start('Linear_Solver_For_Multiphase')
      if(solver .eq. 'ACM') then
        vof % tol   = PICO
        call Acm(sol,            &
                 vof % n,        &
                 b,              &
                 vof % precond,  &
                 vof % niter,    &     ! number of V cycles
                 vof % tol,      &
                 vof % res)
        stop
      else
        call Cg(sol,            &
                vof % n,        &
                b,              &
                vof % precond,  &
                vof % niter,    &      ! max number of iterations
                exec_iter,      &      ! executed number of iterations
                vof % tol,      &
                vof % res)
      end if
      call Cpu_Timer_Mod_Stop('Linear_Solver_For_Multiphase')

      call Info_Mod_Iter_Fill_At(1, 6, vof % name, exec_iter, vof % res)

      !---------------------------!
      !   Correct Beta at faces   !
      !---------------------------!
      call Multiphase_Mod_Vof_Correct_Beta(vof,           &
                                           vol_flux % n,  &
                                           beta_f,        &
                                           dt)

      !------------------------!
      !   Correct Boundaries   !
      !------------------------!
      do s = 1, grid % n_faces

        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)

        if (c2 < 0) then
          if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW) then
            vof % n(c2) = vof % n(c1)
          end if
        end if

      end do

      !---------------------!
      !   Noise Reduction   !
      !---------------------!
      do c = 1, grid % n_cells
        if(abs(vof % n(c)) < TINY) then
          vof % n(c) = 0.0
        end if

        if(abs(vof % n(c) - 1.0) < TINY) then
          vof % n(c) = 1.0
        end if
      end do

      corr_cicsam = .true.
      loop_look: do c = 1, grid % n_cells
        if(vof % n(c) < 0.0) then
          corr_cicsam = .false.
          exit loop_look
        end if
      end do loop_look

      if(corr_cicsam) then
        exit loop_corr
      end if
    end do loop_corr

    !-----------------------------!
    !   Correct Volume Fraction   !
    !-----------------------------!
    do c = 1, grid % n_cells
      if(vof % n(c) < 0.0) then
        vof % n(c) = 0.0
      end if

      if(vof % n(c) > 1.0) then
        vof % n(c) = 1.0
      end if
    end do

  end if


  !-----------------------!
  !   Update properties   !
  !-----------------------!
  do c=1,grid % n_cells
    density(c) = vof % n(c)         * phase_dens(1)         &
               + (1.0 - vof % n(c)) * phase_dens(2)
    viscosity(c) = vof % n(c)         * phase_visc(1)       &
                 + (1.0 - vof % n(c)) * phase_visc(2)
  end do

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

        if (vol_flux % n(s)>=0.0) then
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
        if (flow % flux(s) >= 0) then
          donor = c1
          accept = c2
        else
          donor = c2
          accept = c1
        end if

        vof_f(s) = 0.5 * ((1.0 - beta_f(s)) * ( vof % n(donor)                &
                                              + vof % o(donor) )              &
                               + beta_f(s)  * ( vof % n(accept)               &
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

  !----------------------!
  !   Density at faces   !
  !----------------------!
  do s=1, grid % n_faces
    dens_face(s) = vof_f(s) * phase_dens(1)  &
                 + (1.0 - vof_f(s)) * phase_dens(2)
  end do

  call Cpu_Timer_Mod_Stop('Compute_Multiphase (without solvers)')

  end subroutine


