!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, turb, Vof, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: swarm
  integer, intent(in)      :: n     ! time step
  real,    intent(in)      :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: i, j, k, n_parts_in_buffers
  real    :: x, y, z, dy, dz, my, mz
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NJ = 16, NK = 4
!==============================================================================!

  !----------------------------------------------------!
  !   Initialize particles only in the 1st time step   !
  !----------------------------------------------------!
  if(n .eq. 1) then

    dy = 4.0 / NJ
    dz = 1.0 / NK

    ! Place particles where you want them
    do j = 1, NJ
      do k = 1, NK
        i = (k-1)*NJ + j  ! particle number

        ! Placing particles (only at the 1st time step)
        x = 0.05
        y = dy * 0.5 + (j-1) * dy
        z = dz * 0.5 + (k-1) * dz

        call random_number(my);  my = (my - 0.5) * dy * 0.4
        call random_number(mz);  mz = (mz - 0.5) * dz * 0.4
        swarm % particle(i) % x_n = x
        swarm % particle(i) % y_n = y + my
        swarm % particle(i) % z_n = z + mz

        swarm % particle(i) % x_o = swarm % particle(k) % x_n
        swarm % particle(i) % y_o = swarm % particle(k) % y_n
        swarm % particle(i) % z_o = swarm % particle(k) % z_n

        ! Searching for the closest cell and node to place the moved particle
        call Swarm_Mod_Find_Nearest_Cell(swarm, i, n_parts_in_buffers)
        call Swarm_Mod_Find_Nearest_Node(swarm, i)
      end do
    end do

  end if

  end subroutine
