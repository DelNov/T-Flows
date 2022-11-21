!==============================================================================!
  subroutine Set_Porosity_In_Cells(Porosity, Grid, reg)
!------------------------------------------------------------------------------!
!   Sets porosity in cells                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Porosity_Type) :: Porosity
  type(Grid_Type)      :: Grid
  integer, intent(in)  :: reg       ! porous region rank
!-----------------------------------[Locals]-----------------------------------!
  type(Stl_Type)            :: stl
  character(len=SL)         :: name_out
  integer                   :: fu, n_facets, v, c, f, l
  real                      :: cx, cy, cz, dot_prod
  real, contiguous, pointer :: por_real(:)   ! just for saving
!==============================================================================!

  call Work % Connect_Real_Cell(por_real)

  ! Open the STL file
  call File % Open_For_Reading_Ascii(Porosity % region(reg) % stl_file, fu)

  !-------------------------------------------------------!
  !   Count all the facets and allocate memory for them   !
  !-------------------------------------------------------!
  rewind(fu)
  n_facets = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') n_facets = n_facets + 1
  end do
  print '(a,i9)', ' # Number of facets on the stl: ', n_facets

  allocate(stl % facet(n_facets))

  !------------------------------------------------!
  !   Read all facet normals, vertex coordinates   !
  !   and compute facet centroids along the way    !
  !------------------------------------------------!
  rewind(fu)
  n_facets = 0
  do
    call File % Read_Line(fu)
    if(Line % tokens(1) .eq. 'endsolid') exit
    if(Line % tokens(1) .eq. 'facet') then
      n_facets = n_facets + 1
      read(Line % tokens(3), *) stl % facet(n_facets) % nx
      read(Line % tokens(4), *) stl % facet(n_facets) % ny
      read(Line % tokens(5), *) stl % facet(n_facets) % nz
      call File % Read_Line(fu)                ! 'outer loop'
      stl % facet(n_facets) % xc = 0.0
      stl % facet(n_facets) % yc = 0.0
      stl % facet(n_facets) % zc = 0.0
      do v = 1, 3
        call File % Read_Line(fu)              ! 'vertex 1, 2 and 3'
        read(Line % tokens(2), *) stl % facet(n_facets) % x(v)
        read(Line % tokens(3), *) stl % facet(n_facets) % y(v)
        read(Line % tokens(4), *) stl % facet(n_facets) % z(v)
        stl % facet(n_facets) % xc = stl % facet(n_facets) % xc    &
                                   + stl % facet(n_facets) % x(v)
        stl % facet(n_facets) % yc = stl % facet(n_facets) % yc    &
                                   + stl % facet(n_facets) % y(v)
        stl % facet(n_facets) % zc = stl % facet(n_facets) % zc    &
                                   + stl % facet(n_facets) % z(v)
      end do
      stl % facet(n_facets) % xc = stl % facet(n_facets) % xc * ONE_THIRD
      stl % facet(n_facets) % yc = stl % facet(n_facets) % yc * ONE_THIRD
      stl % facet(n_facets) % zc = stl % facet(n_facets) % zc * ONE_THIRD
      call File % Read_Line(fu)                ! 'endloop'
    end if
  end do
  print '(a)', ' # Read all stl facets!'

  do c = 1, Grid % n_cells

    ! Assume cell is in the stl
    Porosity % region(reg) % cell_porous(c) = .true.

    do f = 1, n_facets

      ! Vector connecting facet centroid with cell centroid
      cx = stl % facet(f) % xc - Grid % xc(c)
      cy = stl % facet(f) % yc - Grid % yc(c)
      cz = stl % facet(f) % zc - Grid % zc(c)

      dot_prod = cx * stl % facet(f) % nx  &
               + cy * stl % facet(f) % ny  &
               + cz * stl % facet(f) % nz

      ! First time this is negative, cell is not in the stl
      if(dot_prod < 0) then
        Porosity % region(reg) % cell_porous(c) = .false.
        goto 1
      end if
    end do
1   continue

  end do

  ! Save for checking
  por_real(1:Grid % n_cells) = 0.0
  do c = 1, Grid % n_cells
    if(Porosity % region(reg) % cell_porous(c)) por_real(c) = 1.0;
  end do
  l = len_trim(Porosity % region(reg) % stl_file)
  name_out = Porosity % region(reg) % stl_file(1:l-4)
  call Grid % Save_Debug_Vtu(name_out,               &
                             scalar_cell=por_real,   &
                             scalar_name='porosity')

  call Work % Disconnect_Real_Cell(por_real)

  end subroutine
