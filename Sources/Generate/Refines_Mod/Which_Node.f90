!==============================================================================!
  integer function Refines_Mod_Which_Node(ref, Grid, c, n)
!------------------------------------------------------------------------------!
!>  Returns the local number (1-8) of node n in cell c.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type)  :: ref   !! type holding information on refinement
  type(Grid_Type)     :: Grid  !! grid being generated (refined here)
  integer, intent(in) :: c     !! cell number
  integer, intent(in) :: n     !! node number
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
!==============================================================================!

  Refines_Mod_Which_Node = 0

  if (c < 0) then
    print *, '# Which node: Cell non existent !'
    return
  end if

  ! Try the node only
  do i=1,8
    if( Grid % cells_n(i,c) .eq. n) then
      goto 1
    end if
  end do

  ! If it failed try his twins also
  do j=1,Grid % twin_n(n,0)
    do i=1,8
      if( Grid % cells_n(i,c) .eq. Grid % twin_n(n,j)) then
        goto 1
      end if
    end do
  end do

  Refines_Mod_Which_Node = 0
  print *, '# Which node: Trouble, node not found !'
  print *, '# x, y, z = ', Grid % xn(n),  &
                           Grid % yn(n),  &
                           Grid % zn(n)
  print *, '# cell    = ', c, ref % cell_level(c)
  return

1 Refines_Mod_Which_Node = i
  return

  end function Refines_Mod_Which_Node
