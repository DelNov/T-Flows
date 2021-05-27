!==============================================================================!
  subroutine Save_Swarm(Results, swarm, time_step, domain)
!------------------------------------------------------------------------------!
!   Writes particles in VTU file format (for VisIt and Paraview)               !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)      :: Results
  type(Swarm_Type), target :: swarm
  integer                  :: time_step
  integer,        optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),     pointer :: grid
  type(Particle_Type), pointer :: part
  integer                      :: k, fu
  character(SL)                :: name_out
  integer                      :: n_remaining_particles
!-----------------------------[Local parameters]-------------------------------!
  character(len= 0)  :: IN_0 = ''           ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(swarm % n_particles < 1) return

  ! Take aliases for the swarm
  grid => swarm % pnt_grid

  !-----------------------------------------!
  !   Only one processor saves the swarm,   !
  !    therefore it has to be refreshed     !
  !-----------------------------------------!
  call Swarm_Mod_Exchange_Particles(swarm)

  !-------------------------------!
  !   Count remaining particles   !
  !-------------------------------!
  n_remaining_particles = 0
  do k = 1, swarm % n_particles
    part => swarm % particle(k)
    if(.not. part % escaped) then
      n_remaining_particles = n_remaining_particles + 1
    end if
  end do

  if(n_remaining_particles .eq. 0) then
    if(this_proc < 2) then
      print *, '# No particles remaining in the domain, nothing to save!'
    end if
    return
  end if

  !----------------------------!
  !                            !
  !   Create .swarm.vtu file   !
  !                            !
  !----------------------------!

  if(this_proc < 2) then

    call File % Set_Name(name_out,             &
                         time_step=time_step,  &
                         appendix ='-swarm',   &
                         extension='.vtu',     &
                         domain=domain)
    call File % Open_For_Writing_Ascii(name_out, fu)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" '//&
                            'byte_order="LittleEndian">'
    write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(fu,'(a,a,i0.0,a)')   &
                IN_2, '<Piece NumberOfPoints="', n_remaining_particles,  &
                           '" NumberOfCells="0">'

    !-----------------------!
    !                       !
    !   Point coordinates   !
    !                       !
    !-----------------------!
    write(fu,'(a,a)') IN_3, '<Points>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                            '="3" format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                  IN_5, part % x_n, part % y_n, part % z_n
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Points>'

    !----------------!
    !                !
    !   Point data   !
    !                !
    !----------------!
    write(fu,'(a,a)') IN_3, '<PointData Scalars="scalars" vectors="velocity">'

    !--------------------!
    !   Particle i.d.s   !
    !--------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Index" ' // &
                            'format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,i9)') IN_5, k
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------!
    !   Closest cells   !
    !-------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="ClosestCell" ' // &
                            'format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,i9)') IN_5, grid % comm % cell_glo(part % cell)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !----------------------------!
    !   Closest boundary cells   !
    !----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" ' // &
                            'Name="ClosestBndCell" format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,i9)') IN_5, grid % comm % cell_glo(part % bnd_cell)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Particle velocities   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" Name="Velocity" ' // &
                            ' NumberOfComponents="3" format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                       &
                   IN_5, part % u, part % v, part % w
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Velocity magnitudes   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" '  //  &
                            ' Name="VelocityMagnitude" '  //  &
                            ' format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                   IN_5, sqrt(part % u**2 + part % v**2 + part % w**2)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------!
    !   Particle diameters   !
    !------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" '  //  &
                            ' Name="ParticleDiameters" '  //  &
                            ' format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)') IN_5, part % d
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Particle processors   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="Processor" ' // &
                            'format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      if(.not. part % escaped) then
        write(fu,'(a,i9)') IN_5, part % proc
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-----------------------!
    !                       !
    !   End of point data   !
    !                       !
    !-----------------------!
    write(fu,'(a,a)') IN_3, '</PointData>'

    !-----------!
    !           !
    !   Cells   !
    !           !
    !-----------!
    write(fu,'(a,a)') IN_3, '<Cells>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                           ' format="ascii">'
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                           ' format="ascii">'
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type="Int64" Name="types"' //  &
                           ' format="ascii">'
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Cells>'

    !------------!
    !            !
    !   Footer   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_2, '</Piece>'
    write(fu,'(a,a)') IN_1, '</UnstructuredGrid>'
    write(fu,'(a,a)') IN_0, '</VTKFile>'
    close(fu)
  end if

  end subroutine
