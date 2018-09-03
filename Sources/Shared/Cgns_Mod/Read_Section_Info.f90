!==============================================================================!
  subroutine Cgns_Mod_Read_Section_Info(base, block, sect)
!------------------------------------------------------------------------------!
!   Read elements connection info for current sect
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: base, block, sect
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id       ! base index number
  integer           :: block_id      ! block index number
  integer           :: sect_id       ! element section index
  character(len=80) :: sect_name     ! name of the Elements_t node
  integer           :: cell_type     ! types of elements in the section
  integer           :: first_cell    ! index of first element
  integer           :: last_cell     ! index of last element
  integer           :: n_bnd         ! index of last boundary element
  integer           :: parent_flag   ! are the parent cells stored (I guess)
  integer           :: error
  integer           :: cnt, bc
!==============================================================================!

  ! Set input parameters
  base_id  = base
  block_id = block
  sect_id  = sect

  ! Get info for an element section
  call Cg_Section_Read_F(file_id,       & !(in )
                         base_id,       & !(in )
                         block_id,      & !(in )
                         sect_id,       & !(in )
                         sect_name,     & !(out)
                         cell_type,     & !(out)
                         first_cell,    & !(out)
                         last_cell,     & !(out)
                         n_bnd,         & !(out)
                         parent_flag,   & !(out)
                         error)           !(out)
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

  ! Consider boundary conditions defined in this block
  do bc = 1, cgns_base(base) % block(block) % n_bnd_conds
    if(index(trim(sect_name), trim(cgns_base(base) % block(block) % bnd_cond(bc) % name), back = .true.) .ne. 0) then

      if(verbose) then
        print *, '#         ---------------------------------'
        print *, '#         Bnd section name:  ', trim(sect_name)
        print *, '#         ---------------------------------'
        print *, '#         Bnd section index: ', sect
        print *, '#         Bnd section type:  ', ElementTypeName(cell_type)
        print *, '#         First cell:        ',  &
                 cgns_base(base) % block(block) % section(sect) % first_cell
        print *, '#         Last cell:         ',  &
                 cgns_base(base) % block(block) % section(sect) % last_cell
      end if

      ! Count boundary cells
      if ( ElementTypeName(cell_type) .eq. 'QUAD_4') cnt_qua = cnt_qua + cnt
      if ( ElementTypeName(cell_type) .eq. 'TRI_3' ) cnt_tri = cnt_tri + cnt

    end if
  end do

  ! Consider only three-dimensional cells / sections
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