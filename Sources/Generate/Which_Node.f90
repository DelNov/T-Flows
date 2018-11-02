!==============================================================================!
  integer function Which_Node(ref, grid, c, n) 
!------------------------------------------------------------------------------!
!   Returns the local number (1-8) of node n in cell c.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod,     only: twin_n
  use Refines_Mod, only: Refines_Type
  use Grid_Mod,    only: Grid_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Refines_Type) :: ref
  type(Grid_Type)    :: grid
  integer            :: c, n
!-----------------------------------[Locals]-----------------------------------!
  integer :: i, j
!==============================================================================!

  Which_Node = 0

  if (c  < 0) then 
    print *, '# Which node: Cell non existent !'
    return
  end if

  ! Try the node only
  do i=1,8
    if( grid % cells_n(i,c) .eq. n) then
      goto 1
    end if
  end do

  ! If it failed try his twins also
  do j=1,twin_n(n,0)
    do i=1,8
      if( grid % cells_n(i,c) .eq. twin_n(n,j)) then
        goto 1
      end if 
    end do
  end do

  Which_Node = 0
  print *, '# Which node: Trouble, node not found !'
  print *, '# x, y, z = ', grid % xn(n),  &
                           grid % yn(n),  &
                           grid % zn(n)
  print *, '# cell    = ', c, ref % cell_level(c)
  return

1 Which_Node = i
  return

  end function Which_Node
