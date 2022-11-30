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
  integer    :: fu, v, f, i, n
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
    n = int4_array(1)

    Stl % n_facets = n
    allocate(Stl % x(3, n));  Stl % x(:,:) = 0.0
    allocate(Stl % y(3, n));  Stl % y(:,:) = 0.0
    allocate(Stl % z(3, n));  Stl % z(:,:) = 0.0
    allocate(Stl % nx  (n));  Stl % nx (:) = 0.0
    allocate(Stl % ny  (n));  Stl % ny (:) = 0.0
    allocate(Stl % nz  (n));  Stl % nz (:) = 0.0

    !------------------------------------------------!
    !   Read all facet normals, vertex coordinates   !
    !------------------------------------------------!
    do f = 1, n
      call File % Read_Binary_Real4_Array(fu, 3)
      Stl % nx(f) = real4_array(1)
      Stl % ny(f) = real4_array(2)
      Stl % nz(f) = real4_array(3)
      do v = 1, 3
        call File % Read_Binary_Real4_Array(fu, 3)
        Stl % x(v, f) = real4_array(1)
        Stl % y(v, f) = real4_array(2)
        Stl % z(v, f) = real4_array(3)
      end do
      read(fu) byte
      read(fu) byte
    end do
    print '(a)', ' # Read all stl facets!'

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
    print '(a,i9)', ' # Number of facets on the stl: ', f

    Stl % n_facets = f
    allocate(Stl % x(3, f));  Stl % x(:,:) = 0.0
    allocate(Stl % y(3, f));  Stl % y(:,:) = 0.0
    allocate(Stl % z(3, f));  Stl % z(:,:) = 0.0
    allocate(Stl % nx  (f));  Stl % nx (:) = 0.0
    allocate(Stl % ny  (f));  Stl % ny (:) = 0.0
    allocate(Stl % nz  (f));  Stl % nz (:) = 0.0

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
        do v = 1, 3
          call File % Read_Line(fu)              ! 'vertex 1, 2 and 3'
          read(Line % tokens(2), *) Stl % x(v, f)
          read(Line % tokens(3), *) Stl % y(v, f)
          read(Line % tokens(4), *) Stl % z(v, f)
        end do
        call File % Read_Line(fu)                ! 'endloop'
      end if
    end do
    print '(a)', ' # Read all stl facets!'

    close(fu)

  end if

  end subroutine

