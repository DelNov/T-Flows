!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Info(base, block, sect)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, sect
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id       ! base index number
  integer           :: block_id      ! block index number
  integer           :: sect_id       ! element section index
  character(len=80) :: sect_name     ! name of the Elements_t node
  character(len=80) :: int_name      ! name of the interface
  integer           :: int_type      ! type of interface.  1-quad, 2-tri, 3-mix
  integer           :: loc_type      ! type of interface.  1-quad, 2-tri, 3-mix
  integer           :: cell_type     ! types of elements in the section
  integer           :: first_cell    ! index of first element
  integer           :: last_cell     ! index of last element
  integer           :: n_bnd         ! index of last boundary element
  integer           :: parent_flag   ! are the parent cells stored (I guess)
  integer           :: error
  integer           :: cnt, bc, int, pos_in_list, i, j, k
  logical           :: found_in_list
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
    print *, '# Failed to read section ', sect, ' info'
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

  !-----------------------------------------------------!
  !   Consider only boundary conditions in this block   !
  !-----------------------------------------------------!

  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds

    ! If point of b.c. is inside this section -> this section is b.c.
    k = 0
    do i = first_cell, last_cell
      do j = 1, cgns_base(base) % block(block) % bnd_cond(bc) % n_nodes

        if( i .eq. &
          cgns_base(base) % block(block) % bnd_cond(bc) % point_list(j) ) then
          cgns_base(base) % block(block) % bnd_cond(bc) % belongs_to_sect(j) = &
            sect_id
          k = k + 1
        end if

      end do
    end do

    if(verbose .and. k .ne. 0) then
      print *, '#         ---------------------------------'
      print *, '#         Bnd section name:  ', trim(sect_name)
      print *, '#         ---------------------------------'
      print *, '#         Bnd section index ', sect
      print *, '#         Bnd section type:  ', ElementTypeName(cell_type)
      print *, '#         Bnd section has:   ', k, ' boundary faces, '
      print *, '#         They belong to ', &
        trim(cgns_base(base) % block(block) % bnd_cond(bc) % name), ' b.c.'
      !print *, '#         First cell:        ',  &
      !  cgns_base(base) % block(block) % section(sect) % first_cell
      !print *, '#         Last cell:         ',  &
      !  cgns_base(base) % block(block) % section(sect) % last_cell
    end if

    ! Count boundary cells
    if ( ElementTypeName(cell_type) .eq. 'QUAD_4') then
      cnt_bnd_qua = cnt_bnd_qua + k
    end if
    if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then
      cnt_bnd_tri = cnt_bnd_tri + k
    end if

  end do

  !--------------------------------------------!
  !   Consider only interfaces in this block   !
  !--------------------------------------------!
  do int = 1, cgns_base(base) % block(block) % n_interfaces
 
    int_name = trim(cgns_base(base) % block(block) % interface(int) % name)
    int_type = cgns_base(base) % block(block) % interface(int) % int_type

    if(index(trim(sect_name), trim(int_name), back = .true.) .ne. 0) then

      ! Count interface cells
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') then
        cnt_int_qua = cnt_int_qua + cnt
        loc_type = 1
      end if
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) then
        cnt_int_tri = cnt_int_tri + cnt
        loc_type = 2
      end if

      ! Add new interface name, if unique
      found_in_list = .false.
      do i = 1, cnt_int
        if (trim(int_name) .eq. trim(interface_names(i))) then
          pos_in_list = i
          found_in_list = .true.
        end if
      end do

      if (.not. found_in_list) then
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
        end if

        ! Retrive unique id of that interface
        cgns_base(base) % block(block) % interface(int) % id = pos_in_list
      end if

      if(verbose) then
        print *, '#         ---------------------------------'
        print *, '#         Interface name:  ', trim(sect_name)
        print *, '#         ---------------------------------'
        print *, '#         Interface section index: ', sect
        print *, '#         Interface type:  ', ElementTypeName(cell_type)
        print *, '#         First cell:        ',  &
          cgns_base(base) % block(block) % section(sect) % first_cell
        print *, '#         Last cell:         ',  &
          cgns_base(base) % block(block) % section(sect) % last_cell
        print *, '#         Marked for deletion:     ', cgns_base(base) % &
          block(block) % interface(int) % marked_for_deletion          
      end if

    end if
  end do

  !------------------------------------------------------!
  !   Consider only three-dimensional cells / sections   !
  !------------------------------------------------------!
  if ( ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) .or.  &
       ( ElementTypeName(cell_type) .eq. 'PENTA_6') .or.  &
       ( ElementTypeName(cell_type) .eq. 'TETRA_4') ) then

    if(verbose) then
      print *, '#         ---------------------------------'
      print *, '#         Cell section name: ', sect_name
      print *, '#         ---------------------------------'
      print *, '#         Cell section idx:  ', sect
      print *, '#         Cell section type: ', ElementTypeName(cell_type)
      print *, '#         First cell:        ',  &
        cgns_base(base) % block(block) % section(sect) % first_cell
      print *, '#         Last cell:         ',  &
        cgns_base(base) % block(block) % section(sect) % last_cell
    end if

    ! Count cells in sect
    if ( ElementTypeName(cell_type) .eq. 'HEXA_8' ) cnt_hex = cnt_hex + cnt
    if ( ElementTypeName(cell_type) .eq. 'PYRA_5' ) cnt_pyr = cnt_pyr + cnt
    if ( ElementTypeName(cell_type) .eq. 'PENTA_6') cnt_wed = cnt_wed + cnt
    if ( ElementTypeName(cell_type) .eq. 'TETRA_4') cnt_tet = cnt_tet + cnt

  end if

  end subroutine