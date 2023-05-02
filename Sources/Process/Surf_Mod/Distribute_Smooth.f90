!==============================================================================!
  subroutine Distribute_Smooth(Surf, smooth)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type),  target :: Surf
  type(Var_Type),    target :: smooth
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Vert_Type),   pointer :: Vert(:)
  integer,           pointer :: nv
  integer                    :: v, c
!==============================================================================!

  ! Take aliases
  Grid => Surf % pnt_grid
  Vert => Surf % Vert
  nv   => Surf % n_verts

  if(Sequential_Run()) then
    do v = 1, nv
      c = Surf % Vert(v) % cell  ! get nearest cell
      Surf % Vert(v) % smooth   = smooth % n(c)
      Surf % Vert(v) % smooth_x = smooth % x(c)
      Surf % Vert(v) % smooth_y = smooth % y(c)
      Surf % Vert(v) % smooth_z = smooth % z(c)
    end do
  else
    do v = 1, nv

      Surf % buff_v(v) = 0.0
      Surf % buff_x(v) = 0.0
      Surf % buff_y(v) = 0.0
      Surf % buff_z(v) = 0.0
      Surf % buff_n(v) = 0

      c = Surf % Vert(v) % cell  ! get nearest cell
      if(c > 0) then             ! if cell is in this processor
        if(Grid % Comm % cell_proc(c) .eq. This_Proc()) then
          Surf % buff_v(v) = smooth % n(c)
          Surf % buff_x(v) = smooth % x(c)
          Surf % buff_y(v) = smooth % y(c)
          Surf % buff_z(v) = smooth % z(c)
          Surf % buff_n(v) = Surf % buff_n(v) + 1
        end if
      end if

    end do

    call Global % Sum_Real_Array(nv, Surf % buff_v)
    call Global % Sum_Real_Array(nv, Surf % buff_x)
    call Global % Sum_Real_Array(nv, Surf % buff_y)
    call Global % Sum_Real_Array(nv, Surf % buff_z)
    call Global % Sum_Int_Array (nv, Surf % buff_n)

    do v = 1, nv
      Surf % Vert(v) % smooth   = Surf % buff_v(v) / Surf % buff_n(v)
      Surf % Vert(v) % smooth_x = Surf % buff_x(v) / Surf % buff_n(v)
      Surf % Vert(v) % smooth_y = Surf % buff_y(v) / Surf % buff_n(v)
      Surf % Vert(v) % smooth_z = Surf % buff_z(v) / Surf % buff_n(v)
    end do
  end if

  end subroutine
