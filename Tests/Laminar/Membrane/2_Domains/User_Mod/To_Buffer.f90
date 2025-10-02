!==============================================================================!
  subroutine Interface_Mod_Exchange(inter, var1, var2, v, boundary)
!------------------------------------------------------------------------------!
!   Copy values of specified variables to buffer.                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type) :: inter
  real                 :: var1(-inter % pnt_grid1 % n_bnd_cells :  &
                                inter % pnt_grid1 % n_cells)
  real                 :: var2(-inter % pnt_grid2 % n_bnd_cells :  &
                                inter % pnt_grid2 % n_cells)
  integer              :: v         ! variable number
  logical, optional    :: boundary  ! exchange boundary values (not inside)
!-----------------------------------[Locals]-----------------------------------!
  integer :: n1, n2, b, c, n, n_tot
  logical :: inside                 ! exchange inside values (not boundary)
!==============================================================================!

  ! Handle optional parameter
  if(present(boundary)) then
    if(      boundary) inside = .false.  ! if boundary true, inside is false
    if(.not. boundary) inside = .true.   ! if boudary false, inside is true
  else
    inside = .true.
  end if

  !----------------------------------------------------!
  !   Copy values from inside cells to the interface   !
  !----------------------------------------------------!

  n_tot = inter % n_tot

  if(n_tot > 0) then

    ! Nullify values at the interface
    inter % phi_1(1:n_tot, v) = 0.0
    inter % phi_2(1:n_tot, v) = 0.0

    ! On the side of domain 1
    do n1 = 1, inter % n1_sub
      c = inter % cell_1(n1)   ! domain 1, cell inside
      b = inter % bcel_1(n1)   ! domain 1, boundary cell
      n = inter % face_1(n1)   ! domain 1, face
      if(inside) then
        inter % phi_1(n, v) = var1(c)
      else
        inter % phi_1(n, v) = var1(b)
      end if
    end do

    ! On the side of domain 2
    do n2 = 1, inter % n2_sub
      c = inter % cell_2(n2)   ! domain 2, cell inside
      b = inter % bcel_2(n2)   ! domain 2, boundary cell
      n = inter % face_2(n2)   ! domain 2, face
      if(inside) then
        inter % phi_2(n, v) = var2(c)
      else
        inter % phi_2(n, v) = var2(b)
      end if
    end do

    ! Here we exchange (global sum) of phi_1 and phi_2
    call Global % Sum_Real_Array(n_tot, inter % phi_1(1:n_tot,v))
    call Global % Sum_Real_Array(n_tot, inter % phi_2(1:n_tot,v))

  end if

  end subroutine
