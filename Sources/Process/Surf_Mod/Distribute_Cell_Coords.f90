!==============================================================================!
  subroutine Distribute_Cell_Coords(Surf)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Vert_Type), pointer :: Vert(:)
  integer,         pointer :: nv
  integer                  :: v, c
!==============================================================================!

  ! Take aliases
  Grid => Surf % pnt_grid
  Vert => Surf % Vert
  nv   => Surf % n_verts

  if(Sequential_Run()) then
    do v = 1, nv
      c = Surf % Vert(v) % cell  ! get nearest cell
      Surf % Vert(v) % cell_x = Grid % xc(c)
      Surf % Vert(v) % cell_y = Grid % yc(c)
      Surf % Vert(v) % cell_z = Grid % zc(c)
    end do
  else
    do v = 1, nv

      Surf % buff_x(v) = 0.0
      Surf % buff_y(v) = 0.0
      Surf % buff_z(v) = 0.0
      Surf % buff_n(v) = 0

      c = Surf % Vert(v) % cell  ! get nearest cell
      if(c > 0) then             ! if cell is in this processor
        if(Cell_In_This_Proc(c)) then
          Surf % buff_x(v) = Grid % xc(c)
          Surf % buff_y(v) = Grid % yc(c)
          Surf % buff_z(v) = Grid % zc(c)
          Surf % buff_n(v) = Surf % buff_n(v) + 1
        end if
      end if

    end do

    call Global % Sum_Real_Array(nv, Surf % buff_x)
    call Global % Sum_Real_Array(nv, Surf % buff_y)
    call Global % Sum_Real_Array(nv, Surf % buff_z)
    call Global % Sum_Int_Array (nv, Surf % buff_n)

    do v = 1, nv
      Surf % Vert(v) % cell_x = Surf % buff_x(v) / Surf % buff_n(v)
      Surf % Vert(v) % cell_y = Surf % buff_y(v) / Surf % buff_n(v)
      Surf % Vert(v) % cell_z = Surf % buff_z(v) / Surf % buff_n(v)
    end do
  end if

  end subroutine
