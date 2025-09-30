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
  integer                   :: c, f, l, n, i_ver, hits
  real                      :: xyz(3), xn(3), yn(3), zn(3), q(3), s(3)
  real                      :: a, g, u, t, mx, my, mz, nx, ny, nz, d(3), h(3)
  real                      :: e1(3), e2(3)
  real                      :: box_x, box_y, box_z
  real                      :: xc, yc, zc, xf, yf, zf, v
  real, contiguous, pointer :: por_real(:)   ! just for saving
!==============================================================================!

  call Work % Connect_Real_Cell(por_real)

  !---------------------------------------------------!
  !                                                   !
  !   An STL file was defined for the porous region   !
  !                                                   !
  !---------------------------------------------------!
  if(Por % region(reg) % stl_name .ne. '') then

    ! Get the bounding box coordinates
    call Grid % Bounding_Box(nx, ny, nz, mx, my, mz)
    box_x = mx - nx
    box_y = my - ny
    box_z = mz - nz

    ! Set the point far far away
    xf = box_x * KILO
    yf = box_y * KILO
    zf = box_z * KILO

    !------------------------------------!
    !   Create an STL object from file   !
    !------------------------------------!
    call Por % region(reg) % Stl % Create_From_File(Por % region(reg) % stl_name)

    !-------------------------------------------!
    !   Use the STL object to define porosity   !
    !- - - - - - - - - - - - - - - - - - - - - -!
    !   Note:                                   !
    !   - T_FLOWS_PROGRAM 1 is Convert          !
    !   - T_FLOWS_PROGRAM 4 is Process          !
    !-------------------------------------------!
# if T_FLOWS_PROGRAM == 4
    do c = Cells_In_Domain_And_Buffers()
# elif T_FLOWS_PROGRAM == 1
    do c = 1, Grid % n_cells
# endif

      ! Cell's center coordinates
      xc = Grid % xc(c)
      yc = Grid % yc(c)
      zc = Grid % zc(c)

      ! Ray/segment direction
      d(1) = xf - xc
      d(2) = yf - yc
      d(3) = zf - zc

      ! Assume cell is in the stl
      Grid % por(c) = 0

      ! Initialize number of hits
      hits = 0

      do f = 1, Por % region(reg) % Stl % N_Facets()

        ! Fetch all vertex coordinates
        do i_ver = 1, 3
          xyz(1:3) = Por % region(reg) % Stl % Facets_Vert_Coords(f, i_ver)
          xn(i_ver) = xyz(1)
          yn(i_ver) = xyz(2)
          zn(i_ver) = xyz(3)
        end do

        e1(1) = xn(2) - xn(1)
        e1(2) = yn(2) - yn(1)
        e1(3) = zn(2) - zn(1)

        e2(1) = xn(3) - xn(1)
        e2(2) = yn(3) - yn(1)
        e2(3) = zn(3) - zn(1)

        h = Math % Cross_Product(d, e2)

        a = dot_product(e1, h)

        ! Check if ray is co-planar with the facet
        if(abs(a) < TINY) cycle

        g = 1.0 / a

        ! S = P - A
        s(1) = xc - xn(1)
        s(2) = yc - yn(1)
        s(3) = zc - zn(1)

        ! u = g * (s · h)
        u = g * dot_product(s, h)

        ! Check if (what on earth are we checking here?)
        if(u .lt. 0.0 .or. u .gt. 1.0) cycle

        ! q = s x e1
        q = Math % Cross_Product(s, e1)

        ! v = g * (D · q)
        v = g * dot_product(d, q)

        ! Check if barycentric v, or u+v
        if(v .lt. 0.0 .or. (u + v) .gt. 1.0) cycle

        ! t = f * (e2 · q)
        t = g * dot_product(e2, q)

        if(t .lt. 0.0) cycle
        if(t .gt. 1.0) cycle

        hits = hits + 1
      end do  ! through facets in the region

      if(mod(hits,2) .eq. 1) then
        Grid % por(c) = reg
        Por % region(reg) % cell_porous(c) = .true.
      end if
    end do    ! through cells

  !-----------------------------------------------------------------------!
  !                                                                       !
  !   Porous regions defined with grid, nothing to be done here, really   !
  !                                                                       !
  !-----------------------------------------------------------------------!
  else
  end if      ! porous region defined in an STL file

  ! Save for checking
  por_real(1:Grid % n_cells) = 0.0
# if T_FLOWS_PROGRAM == 4
  do c = Cells_In_Domain_And_Buffers()
# elif T_FLOWS_PROGRAM == 1
  do c = 1, Grid % n_cells
# endif
    por_real(c) = Grid % por(c)
  end do
  l = len_trim(Por % region(reg) % Stl % name)
  name_out = 'porosity-field'
  call Grid % Save_Debug_Vtu(name_out,               &
                             scalar_cell=por_real,   &
                             scalar_name='porosity')

  call Work % Disconnect_Real_Cell(por_real)

  end subroutine
