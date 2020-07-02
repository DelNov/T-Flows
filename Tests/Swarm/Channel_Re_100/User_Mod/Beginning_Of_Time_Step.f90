!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(flow, turb, mult, swarm, n, time)
!------------------------------------------------------------------------------!
!   This function is called at the end of time step.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: n     ! time step
  real                          :: time  ! physical time
!----------------------------------[Locals]------------------------------------!
  integer :: i, j, k, n_parts_in_buffers
  real    :: x, y, z, dy, dz
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: NJ = 5
  integer, parameter :: NZ = 5
!==============================================================================!

  !-------------------!
  !   1st time step   !
  !-------------------!
  if(n .eq. 1) then

    dy = 1.0 / NJ
    dz = 1.0 / NZ

    ! Initializing both deposition and departure counters
    swarm % cnt_d = 0
    swarm % cnt_e = 0
    swarm % cnt_r = 0

    ! Place 100 particles where you want them
    do j = 1, NJ
      do k = 1, NZ
        i = (k-1)*NJ + j  ! particle number

        ! Placing particles (only at the 1st time step)
        x = 0.05
        y = dy * 0.5 + (j-1) * dy
        z = dz * 0.5 + (k-1) * dz

        swarm % particle(i) % x_n = x
        swarm % particle(i) % y_n = y
        swarm % particle(i) % z_n = z

        ! Searching for the closest cell and node to place the moved particle
        call Swarm_Mod_Find_Nearest_Cell(swarm, i, n_parts_in_buffers)
        call Swarm_Mod_Find_Nearest_Node(swarm, i)
      end do
    end do

  end if

  end subroutine
