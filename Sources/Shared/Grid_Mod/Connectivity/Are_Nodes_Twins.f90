!==============================================================================!
  logical function Are_Nodes_Twins(Grid, n1, n2)
!------------------------------------------------------------------------------!
!   Checks if the nodes are twins, i.e. are they shared on periodicity         !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Grid_Type)     :: Grid    !! parent class
  integer,  intent(in) :: n1, n2  !! node
!-----------------------------------[Locals]-----------------------------------!
  integer :: n, check
!==============================================================================!

  ! This function is used only from Generate
  Assert(PROGRAM_NAME .eq. 'Generate')

  check = 0

  do n = 1, Grid % twin_n(n1,0)
    if(Grid % twin_n(n1,n) .eq. n2) check = check + 1
  end do

  do n = 1, Grid % twin_n(n2,0)
    if(Grid % twin_n(n2,n) .eq. n1) check = check + 1
  end do

  if(check .eq. 2) then
    Are_Nodes_Twins = .true.
    return
  else if(check .eq. 0) then
    Are_Nodes_Twins = .false.
    return
  else
    print *, '# ERROR in Are_Nodes_Twins!  Major trouble !  Stopping !'
    stop
  end if

  end function
