!==============================================================================!
  subroutine Set_Porosity_In_Cells(Por, Grid, reg)
!------------------------------------------------------------------------------!
!   Sets porosity in cells                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Porosity_Type) :: Por
  type(Grid_Type)      :: Grid
  integer, intent(in)  :: reg       ! porous region rank
!-----------------------------------[Locals]-----------------------------------!
  character(len=SL)         :: name_out
  integer                   :: c, f, l
  real                      :: xyz(3), n(3)
  real                      :: nx, ny, nz, xc, yc, zc, cx, cy, cz, dot_prod
  real, contiguous, pointer :: por_real(:)   ! just for saving
!==============================================================================!

  call Work % Connect_Real_Cell(por_real)

  !------------------------------------!
  !   Create an STL object from file   !
  !------------------------------------!
  call Por % region(reg) % Stl % Create_From_File(Por % region(reg) % stl_name)

  !-------------------------------------------!
  !   Use the STL object to define porosity   !
  !-------------------------------------------!
  do c = 1, Grid % n_cells

    ! Assume cell is in the STL
    Por % region(reg) % cell_porous(c) = .true.

    do f = 1, Por % region(reg) % Stl % N_Facets()

      ! Compute facet's center of gravity
      xyz(1:3) = Por % region(reg) % Stl % Facet_Coords(f)
      xc = xyz(1);  yc = xyz(2);  zc = xyz(3)

      n(1:3) = Por % region(reg) % Stl % Facet_Normal(f)
      nx = n(1);  ny = n(2);  nz = n(3)

      ! Vector connecting facet centroid with cell centroid
      cx = xc - Grid % xc(c)
      cy = yc - Grid % yc(c)
      cz = zc - Grid % zc(c)

      dot_prod = cx * nx + cy * ny + cz * nz

      ! First time this is negative, cell is not in the stl
      if(dot_prod < 0) then
        Por % region(reg) % cell_porous(c) = .false.
        goto 1
      end if
    end do
1   continue

  end do

  ! Save for checking
  por_real(1:Grid % n_cells) = 0.0
  do c = 1, Grid % n_cells
    if(Por % region(reg) % cell_porous(c)) por_real(c) = 1.0;
  end do
  l = len_trim(Por % region(reg) % Stl % name)
  name_out = Por % region(reg) % Stl % name(1:l-4)
  call Grid % Save_Debug_Vtu(name_out,               &
                             scalar_cell=por_real,   &
                             scalar_name='porosity')

  call Work % Disconnect_Real_Cell(por_real)

  end subroutine
