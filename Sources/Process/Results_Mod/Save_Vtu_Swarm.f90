!==============================================================================!
  subroutine Save_Vtu_Swarm(Results, Swarm, time_step, domain)
!------------------------------------------------------------------------------!
!   Writes particles in VTU file format (for VisIt and Paraview)               !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)      :: Results
  type(Swarm_Type), target :: Swarm
  integer                  :: time_step
  integer,        optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type),     pointer :: Grid
  type(Field_Type),    pointer :: Flow
  type(Particle_Type), pointer :: Part
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
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Results)
!==============================================================================!

  if(Swarm % n_particles < 1) return

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Take aliases for the Swarm
  Grid => Swarm % pnt_grid
  Flow => Swarm % pnt_flow

  !-----------------------------------------!
  !   Only one processor saves the Swarm,   !
  !    therefore it has to be refreshed     !
  !-----------------------------------------!
  call Swarm_Mod_Exchange_Particles(Swarm)

  !-------------------------------!
  !   Count remaining particles   !
  !-------------------------------!
  n_remaining_particles = 0
  do k = 1, Swarm % n_particles
    Part => Swarm % Particle(k)
    if(.not. Part % escaped) then
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
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' NumberOfComponents="3"'   //  &
                            ' format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                  IN_5, Part % x_n, Part % y_n, Part % z_n
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_3, '</Points>'

    !----------------!
    !                !
    !   Point data   !
    !                !
    !----------------!
    write(fu,'(a,a)') IN_3, '<PointData>'

    !--------------------!
    !   Particle i.d.s   !
    !--------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Index" ' // &
                            'format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,i9)') IN_5, k
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------!
    !   Closest cells   !
    !-------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="ClosestCell" '    //  &
                            'format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,i9)') IN_5, Grid % Comm % cell_glo(Part % cell)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !----------------------------!
    !   Closest boundary cells   !
    !----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="ClosestBndCell" format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,i9)') IN_5, Grid % Comm % cell_glo(Part % bnd_cell)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Particle velocities   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' Name="Velocity" '         // &
                            ' NumberOfComponents="3" format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                   IN_5, Part % u, Part % v, Part % w
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Velocity magnitudes   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp    //  &
                            ' Name="VelocityMagnitude" '  //  &
                            ' format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')  &
                   IN_5, sqrt(Part % u**2 + Part % v**2 + Part % w**2)
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------!
    !   Particle diameters   !
    !------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp    //  &
                            ' Name="ParticleDiameters" '  //  &
                            ' format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)') IN_5, Part % d
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !----------------------!
    !   Particle density   !
    !----------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' Name="ParticleDensity" '  //  &
                            ' format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)') IN_5, Part % density
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------!
    !   Fluid density   !
    !-------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//floatp  //  &
                            ' Name="FluidDensity" '     //  &
                            ' format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)') IN_5, Part % dens_fluid
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Particle processors   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="Processor" '      //  &
                            'format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        write(fu,'(a,i9)') IN_5, Part % proc
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !------------------------------!
    !   Particle deposited state   !
    !------------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="Deposited" '      //  &
                            'format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        if(Part % deposited) then
          write(fu,'(a,i9)') IN_5, 1
        else
          write(fu,'(a,i9)') IN_5, 0
        end if
      end if
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !----------------------------!
    !   Particle trapped state   !
    !----------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="Trapped" '  //  &
                            'format="ascii">'
    do k = 1, Swarm % n_particles
      Part => Swarm % Particle(k)
      if(.not. Part % escaped) then
        if(Part % trapped) then
          write(fu,'(a,i9)') IN_5, 1
        else
          write(fu,'(a,i9)') IN_5, 0
        end if
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
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp  //  &
                            ' Name="connectivity"'    //  &
                           ' format="ascii">'
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="offsets"' //  &
                           ' format="ascii">'
    write(fu,'(a,a)') IN_4, '</DataArray>'
    write(fu,'(a,a)') IN_4, '<DataArray type='//intp//' Name="types"' //  &
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
