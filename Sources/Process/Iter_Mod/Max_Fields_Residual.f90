!==============================================================================!
  pure real function Max_Fields_Residual(Iter, Flow, n_dom)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Iter_Type), intent(in) :: Iter
  type(Field_Type), intent(in) :: Flow(MD)
  integer,          intent(in) :: n_dom
!-----------------------------------[Locals]-----------------------------------!
  integer :: d, sc
  real    :: mr
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

