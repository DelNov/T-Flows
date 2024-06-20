!==============================================================================!
  subroutine Update_Boundary_Values(Process, Flow, Grid, update)
!------------------------------------------------------------------------------!
!>  This is s a simplified version from the same subroutine in Process_Cpu
!>  as it reads only boundary conditions releated to momentum and enthalpy
!>  conservation equations.  Hopefully, as more modules are ported to
!>  Process_Gpu, this source file will get closer and closer to its origin
!>  from Process_Cpu.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)      :: Process  !! parent class
  type(Field_Type), target :: Flow     !! flow object
  type(Grid_Type)          :: Grid
  character(*)             :: update   !! character switch to control
                                       !! which variables to update
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Update_Boundary_Values')

  call String % To_Upper_Case(update)

  if(update .ne. 'MOMENTUM' .and.  &
     update .ne. 'ALL') then
    call Message % Error(72,                                              &
                         'Invalid parameter in function call. Exiting!',  &
                         file=__FILE__, line=__LINE__, one_proc=.true.)
  end if

  !--------------!
  !              !
  !   Momentum   !
  !              !
  !--------------!
  if( (update .eq. 'MOMENTUM' .or. update .eq. 'ALL') ) then

    ! On the boundary perform the extrapolation
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
         Grid % region % type(reg) .eq. PRESSURE .or.  &
         Grid % region % type(reg) .eq. SYMMETRY) then

        !$acc parallel loop                                                  &
        !$acc present(grid_region_f_face, grid_region_l_face, grid_faces_c,  &
        !$acc         flow_u_n, flow_v_n, flow_w_n)
        do s = Faces_In_Region_Gpu(reg)
          c1 = grid_faces_c(1,s)  ! inside cell
          c2 = grid_faces_c(1,s)  ! boundary cell

          flow_u_n(c2) = flow_u_n(c1)
          flow_v_n(c2) = flow_v_n(c1)
          flow_w_n(c2) = flow_w_n(c1)
        end do  ! faces
        !$acc end parallel

      end if    ! boundary condition
    end do      ! region

  end if  ! update momentum

  call Profiler % Stop('Update_Boundary_Values')

  end subroutine
