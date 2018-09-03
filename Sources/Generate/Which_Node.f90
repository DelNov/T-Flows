!==============================================================================!
  integer function Which_Node(grid, c, n) 
!------------------------------------------------------------------------------!
!   Returns the local number (1-8) of node n in cell c.                        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Gen_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: c, n
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
  print *, '# cell    = ', c, level(c)
  return

1 Which_Node = i
  return

  end function Which_Node
