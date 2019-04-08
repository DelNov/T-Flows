!==============================================================================!
  subroutine Save_Swarm(swarm, name_save)
!------------------------------------------------------------------------------!
!   Writes particles in VTU file format (for VisIt and Paraview)               !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Swarm_Type), target :: swarm
  character(len=*)         :: name_save
!----------------------------------[Locals]------------------------------------!
  type(Particle_Type), pointer :: part
  integer                      :: k
  character(len=80)            :: name_out_9, store_name
!-----------------------------[Local parameters]-------------------------------!
  character(len= 0)  :: IN_0 = ''           ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(swarm % n_particles < 1) return

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  !----------------------------!
  !                            !
  !   Create .swarm.vtu file   !
  !                            !
  !----------------------------!

  if(this_proc < 2) then

    call Name_File(0, name_out_9, '.swarm.vtu')

    open(9, file=name_out_9)
    print *, '# Creating file: ', trim(name_out_9)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //&
                           'byte_order="LittleEndian">'
    write(9,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(9,'(a,a,i0.0,a)')   &
               IN_2, '<Piece NumberOfPoints="', swarm % n_particles,  &
                          '" NumberOfCells="0">'

    !-----------------------!
    !                       !
    !   Point coordinates   !
    !                       !
    !-----------------------!
    write(9,'(a,a)') IN_3, '<Points>'
    write(9,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                           '="3" format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      write(9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                 IN_5, part % x_n, part % y_n, part % z_n
    end do
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_3, '</Points>'

    !-----------!
    !           !
    !   Cells   !
    !           !
    !-----------!
    write(9,'(a,a)') IN_3, '<Cells>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                           ' format="ascii">'
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                           ' format="ascii">'
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="types"' //  &
                           ' format="ascii">'
    write(9,'(a,a)') IN_4, '</DataArray>'
    write(9,'(a,a)') IN_3, '</Cells>'

    !------------!
    !            !
    !   Footer   !
    !            !
    !------------!
    write(9,'(a,a)') IN_2, '</Piece>'
    write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
    write(9,'(a,a)') IN_0, '</VTKFile>'
    close(9)
  end if

  ! Restore the name
  problem_name = store_name

  end subroutine
