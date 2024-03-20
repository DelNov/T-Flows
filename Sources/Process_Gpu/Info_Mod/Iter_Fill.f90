!==============================================================================!
  subroutine Iter_Fill(Info, n)
!------------------------------------------------------------------------------!
!>  The Iter_Fill subroutine is responsible for populating the iteration
!>  information box with basic data, such as the current iteration number and a
!>  placeholder for the mass residual. This function sets up the primary
!>  structure of the iteration info box, which is later populated with more
!>  detailed information about iterative outcomes and targets of various fields
!>  in the simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Info_Type)    :: Info  !! parent, singleton object Info
  integer, intent(in) :: n     !! number of iterations performed in linear solver
!-----------------------------------[Locals]-----------------------------------!
  integer :: c
!==============================================================================!

  if(First_Proc()) then

    ! Write basic info
    c = 1  ! a column
    write(Info % iter % line_iter((c-1)*L_BOX+58 :  &
                                  (c-1)*L_BOX+67),  '(a10)') 'Iteration:'
    write(Info % iter % line_iter((c-1)*L_BOX+68 :  &
                                  (c-1)*L_BOX+70),   '(i3)') n

    c = 5  ! a column
    write(Info % iter % line(1)((c-1)*L_BOX+3 :  &
                                (c-1)*L_BOX+7),    '(a5)') ' VOL:'

  end if

  end subroutine
