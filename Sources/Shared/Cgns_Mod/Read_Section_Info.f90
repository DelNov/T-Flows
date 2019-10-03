!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Info(base, block, sect)
!------------------------------------------------------------------------------!
!   Read sect info, assign b.c, interfaces and count faces, cell               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, sect
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id       ! base index number
  integer           :: block_id      ! block index number
  integer           :: sect_id       ! element section index
  character(len=80) :: sect_name     ! name of the elements_t node
  character(len=80) :: int_name      ! name of the interface
  character(len=80) :: bnd_name      ! name of the boundary section
  integer           :: min_name_l
  integer           :: int_type      ! type of interface.  1-quad, 2-tri, 3-mix
  integer           :: loc_type      ! type of interface.  1-quad, 2-tri, 3-mix
  integer           :: cell_type     ! types of elements in the section
  integer           :: first_cell    ! index of first element
  integer           :: last_cell     ! index of last element
  integer           :: n_bnd         ! index of last boundary element
  integer           :: parent_flag   ! are the parent cells stored (I guess)
  integer           :: error
  integer           :: cnt, bc, int, i, j, f
  integer           :: bc_found, int_found, pos_in_list
  logical           :: found_in_global_int_list
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Get info for an element section
  call Cg_Section_Read_F(file_id,      & !(in )
                         base_id,      & !(in )
                         block_id,     & !(in )
                         sect_id,      & !(in )
                         sect_name,    & !(out)
                         cell_type,    & !(out)
                         first_cell,   & !(out)
                         last_cell,    & !(out)
                         n_bnd,        & !(out)
                         parent_flag,  & !(out)
                         error)          !(out)
  if (error.ne.0) then
    print "(a)", " # Failed to read section ", sect, " info"
    call Cg_Error_Exit_F()
  endif

  ! Fetch received parameters
  cgns_base(base) % block(block) % section(sect) % name        = trim(sect_name)
  cgns_base(base) % block(block) % section(sect) % cell_type   = cell_type
  cgns_base(base) % block(block) % section(sect) % first_cell  = first_cell
  cgns_base(base) % block(block) % section(sect) % last_cell   = last_cell
  cgns_base(base) % block(block) % section(sect) % parent_flag = parent_flag

  ! Number of cells in this section
  cnt = last_cell - first_cell + 1 ! cells in this sections

  !---------------------------------------!
  !   Consider only 2d-cells / sections   !
  !---------------------------------------!

  if ( ElementTypeName(cell_type) .eq. 'QUAD_4' .or. &
       ElementTypeName(cell_type) .eq. 'TRI_3' ) then

    ! For each face
    do i = first_cell, last_cell

      !-----------------------------------------------------!
      !   Assign points from boundary conditions to sect    !
      !-----------------------------------------------------!

      ! Search through all b.c. conditions
      do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

        do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

          ! Take point from b.c. list
          f = cgns_base(base) % block(block) % bnd_cond(bc) % point_list(j)

          if ( i.eq.f ) then
            cgns_base(base) % block(block) % bnd_cond(bc) % &
              belongs_to_sect(j) = sect_id
          end if ! i.eq.f
        end do ! j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
      end do ! bc = 1, cgns_base(base) % block(block) % n_bnd_conds

      !-------------------------------------------!
      !   Assign points from interfaces to sect   !
      !-------------------------------------------!

      ! Search through all interfaces
      do int = 1, cgns_base(base) % block(block) % n_interfaces

        do j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes

          ! Take point from interface list
          f = cgns_base(base) % block(block) % interface(int) % point_list(j)

          if ( i.eq.f ) then
            cgns_base(base) % block(block) % interface(int) % &
              belongs_to_sect(j) = sect_id
          end if ! i.eq.f
        end do ! j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes
      end do !int = 1, cgns_base(base) % block(block) % n_interfaces

    end do ! i = first_cell, last_cell

    !----------------------!
    !   Count b.c. faces   !
    !----------------------!
    do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

      bc_found = 0
      do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes
        if (cgns_base(base) % block(block) % bnd_cond(bc) % &
          belongs_to_sect(j) .eq. sect_id) then
          bc_found = bc_found + 1
        end if
      end do ! j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

      if(bc_found .ne. 0) then

        ! Count boundary cells
        if ( ElementTypeName(cell_type) .eq. 'QUAD_4') then
          cnt_bnd_qua = cnt_bnd_qua + bc_found
        end if
        if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then
          cnt_bnd_tri = cnt_bnd_tri + bc_found
        end if

        if(verbose) then
          print "(a)", " #====================================================="
          print "(a,a28)", " #       Bnd section name: ", trim(sect_name)
          print "(a)", " #-----------------------------------------------------"
          print "(a,i28)", " #       Bnd section index ", sect
          print "(a,a28)", &
            " #       Bnd section type: ", trim(ElementTypeName(cell_type))
          print "(a,i21)", " #       Bnd section has # faces: ", bc_found
          print "(a,a25)", " #       They belong to b.c.: ", &
            trim(cgns_base(base) % block(block) % bnd_cond(bc) % name)
          print "(a,i34)", " #       First cell: ",  &
            cgns_base(base) % block(block) % section(sect) % first_cell
          print "(a,i35)", " #       Last cell: ",  &
            cgns_base(base) % block(block) % section(sect) % last_cell
        end if ! verbose
      end if !bc_found .ne. 0
    end do ! bc = 1, cgns_base(base) % block(block) % n_bnd_conds

    !-------------------------------------------------------!
    !   Count interface faces and mark a few for deletion   !
    !-------------------------------------------------------!
    do int = 1, cgns_base(base) % block(block) % n_interfaces

      int_found = 0
      do j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes
        if (cgns_base(base) % block(block) % interface(int) % &
          belongs_to_sect(j) .eq. sect_id) then
          int_found = int_found + 1
        end if ! if face is in sect
      end do ! j = 1, cgns_base(base) % block(block) % interface(int) % n_nodes

      if (int_found .ne. 0) then

        ! Count interface cells
        if ( ElementTypeName(cell_type) .eq. 'QUAD_4') then
          cnt_int_qua = cnt_int_qua + int_found
          loc_type = 1
        end if
        if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then
          cnt_int_tri = cnt_int_tri + int_found
          loc_type = 2
        end if

        int_name = trim(cgns_base(base) % block(block) % interface(int) % name)
        int_type = cgns_base(base) % block(block) % interface(int) % int_type

        ! Add new interface name, if unique
        found_in_global_int_list = .false.
        do i = 1, cnt_int
          if (trim(int_name) .eq. trim(interface_names(i))) then
            pos_in_list = i
            found_in_global_int_list = .true.
          end if
        end do ! i = 1, cnt_int

        if (.not. found_in_global_int_list) then
          ! Increase number of unique interfaces
          cnt_int = cnt_int + 1
          interface_names(cnt_int) = trim(int_name)
          cgns_base(base) % block(block) % interface(int) % id = cnt_int
          cgns_base(base) % block(block) % interface(int) % int_type = loc_type
        else
          ! If this interfaces in a list, but defined with different type
          if (loc_type .ne. int_type .and. int_type > 0 .and. int_type < 3) then
            cgns_base(base) % block(block) % interface(int) % int_type = 3 ! mix
          else
            ! This interface name was already added, mark for deletion
            cgns_base(base) % block(block) %  &
              interface(int) % marked_for_deletion = .true.
          end if ! type

          ! Retrive unique id of that interface
          cgns_base(base) % block(block) % interface(int) % id = pos_in_list
        end if ! .not. found_in_global_int_list

        if(verbose) then
         print "(a)",  " #-----------------------------------------------------"
          print "(a,a22)", " #       Interface section name: ", trim(sect_name)
         print "(a)",  " #-----------------------------------------------------"
          print "(a,i21)", " #       Interface section index: ", sect
          print "(a,a29)", &
            " #       Interface type:  ", trim(ElementTypeName(cell_type))
          print "(a,i34)", " #       First cell: ",  &
            cgns_base(base) % block(block) % section(sect) % first_cell
          print "(a,i35)", " #       Last cell: ",  &
            cgns_base(base) % block(block) % section(sect) % last_cell
          print "(a,l25)", " #       Marked for deletion: ", cgns_base(base) % &
            block(block) % interface(int) % marked_for_deletion
        end if ! verbose
      end if ! int_found .ne. 0
    end do ! int = 1, cgns_base(base) % block(block) % n_interfaces

  end if ! QUAD_4, TRI_3

  !------------------------------------------------------!
  !   Consider only three-dimensional cells / sections   !
  !------------------------------------------------------!
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    if(verbose) then
      print "(a)",     " #-----------------------------------------------------"
      print "(a,a24)", " #       3d cell section name: ", trim(sect_name)
      print "(a)",     " #-----------------------------------------------------"
      print "(a,i28)", " #       Cell section idx: ", sect
      print "(a,a27)", " #       Cell section type: ", &
        trim(ElementTypeName(cell_type))
      print "(a,i34)", " #       First cell: ",  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print "(a,i35)", " #       Last cell: ",  &
        cgns_base(base) % block(block) % section(sect) % last_cell
    end if

    ! Globalized first and last cell of this 3d section
    if (pos_of_last_3d_cell .eq. 0) then
      cgns_base(base) % block(block) % section(sect) % first_cell = 1
      cgns_base(base) % block(block) % section(sect) % last_cell  = cnt
      pos_of_last_3d_cell = cnt
    else
      cgns_base(base) % block(block) % section(sect) % first_cell =  &
        pos_of_last_3d_cell + 1
      cgns_base(base) % block(block) % section(sect) % last_cell  =  &
        pos_of_last_3d_cell + cnt
      pos_of_last_3d_cell = pos_of_last_3d_cell + cnt
    end if ! pos_of_last_3d_cell .eq. 0

    if(verbose) then
      print "(a,i24)", " #       Corrected first cell: ",  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print "(a,i25)", " #       Corrected last cell: ",  &
        cgns_base(base) % block(block) % section(sect) % last_cell
    end if ! verbose

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

  end if ! HEXA_8, PYRA_5, PENTA_6, TETRA_4

  end subroutine