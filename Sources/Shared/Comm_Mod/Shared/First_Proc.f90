!==============================================================================!
  pure logical function First_Proc()
!------------------------------------------------------------------------------!
!>  Simple but essential utility in parallel computing environments. It checks
!>  whether the current processor is the first in the sequence. In parallel
!>  executions, it identifies if the processor has a rank of 1, indicating it's
!>  the first in the order. In sequential runs, however, where only one
!>  processor is involved, it checks if the processor's rank is 0.
!>  This functionality is crucial for tasks that need to be executed only from
!>  one processor, quite often limiting the printed output.
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  First_Proc = Global % this_processor < 2

  end function
