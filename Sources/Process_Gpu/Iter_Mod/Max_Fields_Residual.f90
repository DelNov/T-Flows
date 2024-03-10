!==============================================================================!
  pure real function Max_Fields_Residual(Iter, Flow, n_dom)
!------------------------------------------------------------------------------!
!>  Calculates and returns the maximum residual across all fields involved
!>  in a simulation.  It is crucial in deciding when to stop the outer
!>  iterations in the SIMPLE or PISO algorithm.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter      !! parent, singleton object Iter
  type(Field_Type), intent(in) :: Flow(MD)  !! all flow domains in simulation
  integer,          intent(in) :: n_dom     !! actual number of domains
!-----------------------------------[Locals]-----------------------------------!
  integer :: d, sc
  real    :: mr
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Iter)
!==============================================================================!

  mr = 0.0

  do d = 1, n_dom
    mr = max(mr, Flow(d) % vol_res)
    mr = max(mr, Flow(d) % u % res)
    mr = max(mr, Flow(d) % v % res)
    mr = max(mr, Flow(d) % w % res)
    if(Flow(d) % heat_transfer) then
      mr = max(mr, Flow(d) % t % res)
    end if
    do sc = 1, Flow(d) % n_scalars
      mr = max(mr, Flow(d) % scalar(sc) % res)
    end do
  end do

  Max_Fields_Residual = mr

  end function

