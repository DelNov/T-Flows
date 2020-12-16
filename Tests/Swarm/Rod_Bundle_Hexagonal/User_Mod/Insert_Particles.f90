!==============================================================================!
  subroutine User_Mod_Insert_Particles(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer, intent(in)           :: n     ! time step
  real,    intent(in)           :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: i, j, k, n_parts_in_buffers
  real    :: x, y, z, xo, yo, zo, dy, dz, my, mz
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NJ = 16, NK = 16
  real,    parameter :: LY = 0.4, LZ = 0.4
!==============================================================================!

  !----------------------------------------------------!
  !   Initialize particles only in the 1st time step   !
  !----------------------------------------------------!
  if(n .eq. 1) then

    ! Leave 10% margin
    xo = -1.0
    yo =  0.1 * LY - LY/2.0
    zo =  0.1 * LZ
    dy =  0.8 * LY / real(NJ-1)
    dz =  0.8 * LZ / real(NK-1)

    ! Place particles where you want them
    do j = 1, NJ
      do k = 1, NK
        i = (k-1)*NJ + j  ! particle number

        ! Placing particles (only at the 1st time step)
        x = xo
        y = yo + (j-1) * dy
        z = zo + (k-1) * dz
        ! y = yo + dy * LY + (j-1) * dy
        ! z = zo + dz * LZ + (k-1) * dz

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
