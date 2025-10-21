!==============================================================================!
  subroutine Update_Boundary_Values(Process, Grid, Flow, Turb, update)
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
  type(Grid_Type),  target :: Grid     !! grid object
  type(Field_Type), target :: Flow     !! flow object
  type(Turb_Type),  target :: Turb     !! turbulence object
  character(*)             :: update   !! character switch to control
                                       !! which variables to update
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Update_Boundary_Values')

  call String % To_Upper_Case(update)

  if(update .ne. 'MOMENTUM'   .and.  &
     update .ne. 'TURBULENCE' .and.  &
     update .ne. 'ENERGY'     .and.  &
     update .ne. 'SCALARS'    .and.  &
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
        !$tf-acc loop begin
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          Flow % u % n(c2) = Flow % u % n(c1)
          Flow % v % n(c2) = Flow % v % n(c1)
          Flow % w % n(c2) = Flow % w % n(c1)
        end do
        !$tf-acc loop end
      end if    ! boundary condition
    end do      ! region

  end if  ! update momentum

  !-------------------!
  !                   !
  !   Heat transfer   !  =--> the new way
  !                   !
  !-------------------!
  if( (update .eq. 'ENERGY' .or. update .eq. 'ALL') .and.  &
      Flow % heat_transfer) then

    ! On the boundary perform the extrapolation
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW  .or.  &
         Grid % region % type(reg) .eq. PRESSURE .or.  &
         Grid % region % type(reg) .eq. SYMMETRY) then
        !$tf-acc loop begin
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)

          Flow % t % n(c2) = Flow % t % n(c1)
        end do
        !$tf-acc loop end
      end if    ! boundary condition
    end do      ! region

  end if  ! update energy and heat transfer

  call Profiler % Stop('Update_Boundary_Values')

  end subroutine
