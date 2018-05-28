!==============================================================================!
  subroutine Write_Link_To_Mesh_In_File(input_name, base, block)
!------------------------------------------------------------------------------!
!   Writes links for GridCoordinates_t and Elements_t nodes in DB              !
!   This allowes to write a mesh once and then use links to it further         !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod, only: this_proc
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(len=*) :: input_name
  integer          :: base, block
!-----------------------------------[Locals]-----------------------------------!
  integer           :: base_id        ! base index number
  integer           :: block_id       ! block index number
  integer           :: sect_id
  character(len=80) :: file     ! name of the FlowSolution_t node
  character(len=80) :: base_name          ! name of the block
  character(len=80) :: block_name          ! name of the block
  character(len=80) :: name_in_file
  character(len=80) :: nodename
  integer           :: cnt
  integer           :: error
!==============================================================================!

  ! Set input parameters
  file       = input_name

  base_id    = base
  block_id   = block
  base_name  = cgns_base(base_id) % name
  block_name = cgns_base(base_id) % block(block_id) % name

  !-------------------------!
  !   Go to Zone_t DB node  !
  !-------------------------!

  call Cg_Goto_F(file_id,   & !(in )
                 base_id,   & !(in )
                 error,     & !(out)
                 "Zone_t",  & !(in )
                 block_id,  & !(in )
                 "end")

  if (error .ne. 0) then
    print *, "# Failed to navigate to: ", "Zone_t"
    call Cg_Error_Exit_F()
  endif

  !-----------------------!
  !                       !
  !   Coordinates block   !
  !                       !
  !-----------------------!

  nodename = "GridCoordinates"
  name_in_file = "/"//trim(base_name)//"/"&
    //trim(block_name)//"/"//trim(nodename)

  ! Create a link in 'file' to 'name_in_file' and name it 'nodename'
  call Cg_Link_Write_F(trim(nodename),      & ! (in )
                       trim(file),          & ! (in )
                       trim(name_in_file),  & ! (in )
                       error)                 ! (out)

  if (error .ne. 0) then
    print *, "# Failed to create a link to ", trim(file)
    call Cg_Error_Exit_F()
  endif

  ! Print some info
  if(verbose .and. this_proc.lt.2) then
    print *, "#     GridCoordinates linked to : ",  &
      trim(file), ":", trim(name_in_file)
  end if

  !-----------------------------!
  !                             !
  !   Cells connections block   !
  !                             !
  !-----------------------------!
  do sect_id = 1, cgns_base(base) % block(block) % n_sects
    ! cells of sect_id
    if (sect_id.eq.1) cnt = cnt_hex
    if (sect_id.eq.2) cnt = cnt_wed
    if (sect_id.eq.3) cnt = cnt_pyr
    if (sect_id.eq.4) cnt = cnt_tet

    ! if section is not empty
    if (cnt.ne.0) then
      nodename = cgns_base(base_id)%block(block_id)%section(sect_id)%name
      name_in_file = "/"//trim(base_name)//"/" &
        //trim(block_name)//"/"//trim(nodename)

      ! Create a link in 'file' to 'name_in_file' and name it 'nodename'
      call Cg_Link_Write_F(trim(nodename),      & ! (in )
                           trim(file),          & ! (in )
                           trim(name_in_file),  & ! (in )
                           error)                 ! (out)

      if (error .ne. 0) then
        print *, "# Failed to create a link to ", trim(file)
        call Cg_Error_Exit_F()
      endif

      ! Print some info
      if(verbose .and. this_proc.lt.2) then
        print *, "#     ", trim(name_in_file), " linked to : ",  &
          trim(file), ":", trim(name_in_file)
      end if

    end if

  end do

  end subroutine