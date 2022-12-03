!==============================================================================!
  subroutine Create_From_File(Stl, file_name)
!------------------------------------------------------------------------------!
!   Creates an Stl object from a file                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Stl_Type) :: Stl
  character(*)    :: file_name
!-----------------------------------[Locals]-----------------------------------!
  integer    :: fu, v, f, i_ver, n
  integer(1) :: byte
!==============================================================================!

  ! Set the STL's name and open the STL file
  Stl % name = file_name

  !----------------------------------!
  !                                  !
  !   File is in the binary format   !
  !                                  !
  !----------------------------------!
  if(File % Is_In_Binary_Format(Stl % name)) then
    call File % Open_For_Reading_Binary(Stl % name, fu)

    ! First 80 characters are reserved and to be skiped
    do n = 1, 80
      read(fu) byte
    end do

    call File % Read_Binary_Int4_Array(fu, 1)
    f = int4_array(1)

    Stl % n_bnd_cells = -1
    Stl % n_cells     = f
    Stl % n_nodes     = Stl % n_cells * 3
    Stl % n_faces     = 0
    Stl % n_edges     = 0
    allocate(Stl % xn(Stl % n_nodes));  Stl % xn(:) = 0.0
    allocate(Stl % yn(Stl % n_nodes));  Stl % yn(:) = 0.0
    allocate(Stl % zn(Stl % n_nodes));  Stl % zn(:) = 0.0
    allocate(Stl % nx(Stl % n_cells));  Stl % nx(:) = 0.0
    allocate(Stl % ny(Stl % n_cells));  Stl % ny(:) = 0.0
    allocate(Stl % nz(Stl % n_cells));  Stl % nz(:) = 0.0
    allocate(Stl % cells_n_nodes(Stl % n_cells));  Stl % cells_n_nodes(:) = 3
    allocate(Stl % cells_n   (3, Stl % n_cells));  Stl % cells_n(:,:)     = 0
    allocate(Stl % new_n        (Stl % n_nodes));  Stl % new_n(:)         = 0
    allocate(Stl % old_n        (Stl % n_nodes));  Stl % old_n(:)         = 0

    !------------------------------------------------!
    !   Read all facet normals, vertex coordinates   !
    !------------------------------------------------!
    do f = 1, Stl % n_cells
      call File % Read_Binary_Real4_Array(fu, 3)
      Stl % nx(f) = real4_array(1)
      Stl % ny(f) = real4_array(2)
      Stl % nz(f) = real4_array(3)
      do i_ver = 1, 3
        v = (f-1) * 3 + i_ver        ! find global vertex number ...
        Stl % cells_n(i_ver, f) = v  ! ... and store it in cells_n
        call File % Read_Binary_Real4_Array(fu, 3)
        Stl % xn(v) = real4_array(1)
        Stl % yn(v) = real4_array(2)
        Stl % zn(v) = real4_array(3)
      end do
      read(fu) byte
      read(fu) byte
    end do
    print '(a)', ' # Read all STL facets!'

  !---------------------------------!
  !                                 !
  !   File is in the ASCII format   !
  !                                 !
  !---------------------------------!
  else
    call File % Open_For_Reading_Ascii(Stl % name, fu)

    !-------------------------------------------------------!
    !   Count all the facets and allocate memory for them   !
    !-------------------------------------------------------!
    rewind(fu)
    f = 0
    do
      call File % Read_Line(fu)
      if(Line % tokens(1) .eq. 'endsolid') exit
      if(Line % tokens(1) .eq. 'facet') f = f + 1
    end do
    print '(a,i9)', ' # Number of facets on the STL: ', f

    Stl % n_bnd_cells = -1
    Stl % n_cells     = f
    Stl % n_nodes     = Stl % n_cells * 3
    Stl % n_faces     = 0
    Stl % n_edges     = 0
    allocate(Stl % xn(Stl % n_nodes));  Stl % xn(:) = 0.0
    allocate(Stl % yn(Stl % n_nodes));  Stl % yn(:) = 0.0
    allocate(Stl % zn(Stl % n_nodes));  Stl % zn(:) = 0.0
    allocate(Stl % nx(Stl % n_cells));  Stl % nx(:) = 0.0
    allocate(Stl % ny(Stl % n_cells));  Stl % ny(:) = 0.0
    allocate(Stl % nz(Stl % n_cells));  Stl % nz(:) = 0.0
    allocate(Stl % cells_n_nodes(Stl % n_cells));  Stl % cells_n_nodes(:) = 3
    allocate(Stl % cells_n   (3, Stl % n_cells));  Stl % cells_n(:,:)     = 0
    allocate(Stl % new_n        (Stl % n_nodes));  Stl % new_n(:)         = 0
    allocate(Stl % old_n        (Stl % n_nodes));  Stl % old_n(:)         = 0

    !------------------------------------------------!
    !   Read all facet normals, vertex coordinates   !
    !   and compute facet centroids along the way    !
    !------------------------------------------------!
    rewind(fu)
    f = 0
    do
      call File % Read_Line(fu)
      if(Line % tokens(1) .eq. 'endsolid') exit
      if(Line % tokens(1) .eq. 'facet') then
        f = f + 1
        read(Line % tokens(3), *) Stl % nx(f)
        read(Line % tokens(4), *) Stl % ny(f)
        read(Line % tokens(5), *) Stl % nz(f)
        call File % Read_Line(fu)                ! 'outer loop'
        do i_ver = 1, 3
          v = (f-1) * 3 + i_ver        ! find global vertex number ...
          Stl % cells_n(i_ver, f) = v  ! ... and store it in facet_v
          call File % Read_Line(fu)              ! 'vertex 1, 2 and 3'
          read(Line % tokens(2), *) Stl % xn(v)
          read(Line % tokens(3), *) Stl % yn(v)
          read(Line % tokens(4), *) Stl % zn(v)
        end do
        call File % Read_Line(fu)                ! 'endloop'
      end if
    end do
    print '(a)', ' # Read all STL facets!'

    close(fu)

  end if

  !------------------------------------------!
  !                                          !
  !   This is brave: merge duplicate nodes   !
  !                                          !
  !------------------------------------------!
  call Stl % Merge_Duplicate_Nodes()

  end subroutine

