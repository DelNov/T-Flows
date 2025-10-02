!==============================================================================!
  subroutine User_Mod_Insert_Particles(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
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
  if(Time % Curr_Dt() .eq. 1) then

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

        call Swarm % Particle(i) % Insert_At(x, y+my, z+mz,       &
                                             n_parts_in_buffers,  &
                                             Flow=Flow)
      end do
    end do

    !--------------------------------------------------!
    !   Update number of particles in the simulation   !
    !--------------------------------------------------!
    Swarm % n_particles = i

  end if

  end subroutine
