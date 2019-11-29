!==============================================================================!
  subroutine Save_Swarm(swarm, plot_inside)
!------------------------------------------------------------------------------!
!   Writes particles in VTU file format (for VisIt and Paraview)               !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Work_Mod, only: v2_calc   => r_cell_01,  &
                      uu_save   => r_cell_02,  &
                      vv_save   => r_cell_03,  &
                      ww_save   => r_cell_04,  &
                      uv_save   => r_cell_05,  &
                      uw_save   => r_cell_06,  &
                      vw_save   => r_cell_07,  &
                      t2_save   => r_cell_08,  &
                      ut_save   => r_cell_09,  &
                      vt_save   => r_cell_10,  &
                      wt_save   => r_cell_11,  &
                      kin_vis_t => r_cell_12,  &
                      phi_save  => r_cell_13
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Swarm_Type),      target  :: swarm
  type(Grid_Type),       pointer :: grid
  type(Field_Type),      target  :: flow
!  integer                  :: time_step
  logical                  :: plot_inside  ! plot results inside?
!----------------------------------[Locals]------------------------------------!
  type(Particle_Type), pointer :: part
  integer                      :: k, fu, f8, f9, c
  character(len=80)            :: name_out
!-----------------------------[Local parameters]-------------------------------!
  character(len= 0)  :: IN_0 = ''           ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  if(swarm % n_particles < 1) return

  ! Take aliases
  grid => flow % pnt_grid

  !----------------------------!
  !                            !
  !   Create .swarm.vtu file   !
  !                            !
  !----------------------------!

  if(this_proc < 2) then

    call File_Mod_Set_Name(name_out, extension='.swarm.vtu')
    call File_Mod_Open_File_For_Writing(name_out, fu)

    !------------!
    !            !
    !   Header   !
    !            !
    !------------!
    write(fu,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(fu,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //&
                            'byte_order="LittleEndian">'
    write(fu,'(a,a)') IN_1, '<UnstructuredGrid>'

    write(fu,'(a,a,i0.0,a)')   &
                IN_2, '<Piece NumberOfPoints="', swarm % n_particles,  &
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
      write(fu, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                  IN_5, part % x_n, part % y_n, part % z_n
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
      write(fu,'(a,i9)') IN_5, k
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    !-------------------------!
    !   Particle velocities   !
    !-------------------------!
    write(fu,'(a,a)') IN_4, '<DataArray type="Float64" Name="Velocity" ' // &
                            ' NumberOfComponents="3" format="ascii">'
    do k = 1, swarm % n_particles
      part => swarm % particle(k)
      write(fu,'(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                         &
                 IN_5, part % u, part % v, part % w
    end do
    write(fu,'(a,a)') IN_4, '</DataArray>'

    write(fu,'(a,a)') IN_3, '</PointData>'


    !-------------------------!
    !   Particle statistics   !
    !-------------------------!

    ! Statistics for large-scale simulations of turbulence
    if(swarm % swarm_statistics) then
      call Save_Vector(grid, IN_4, IN_5, "MeanParticleVelocity", plot_inside,   &
                                         swarm % u_mean(-grid % n_bnd_cells),   &
                                         swarm % v_mean(-grid % n_bnd_cells),   &
                                         swarm % w_mean(-grid % n_bnd_cells),   &
                                         f8, f9)
      uu_save(:) = 0.0
      vv_save(:) = 0.0
      ww_save(:) = 0.0
      uv_save(:) = 0.0
      uw_save(:) = 0.0
      vw_save(:) = 0.0
      do c = 1, grid % n_cells
        uu_save(c) = swarm % uu(c) - swarm % u_mean(c) * swarm % u_mean(c)
        vv_save(c) = swarm % vv(c) - swarm % v_mean(c) * swarm % v_mean(c)
        ww_save(c) = swarm % ww(c) - swarm % w_mean(c) * swarm % w_mean(c)
        uv_save(c) = swarm % uv(c) - swarm % u_mean(c) * swarm % v_mean(c)
        uw_save(c) = swarm % uw(c) - swarm % u_mean(c) * swarm % w_mean(c)
        vw_save(c) = swarm % vw(c) - swarm % v_mean(c) * swarm % w_mean(c)
      end do
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressXX", plot_inside,  &
                                         uu_save(-grid % n_bnd_cells),         &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressYY", plot_inside,  &
                                         vv_save(-grid % n_bnd_cells),         &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressZZ", plot_inside,  &
                                         ww_save(-grid % n_bnd_cells),         &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressXY", plot_inside,  &
                                         uv_save(-grid % n_bnd_cells),         &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressXZ", plot_inside,  &
                                         uw_save(-grid % n_bnd_cells),         &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanParticleStressYZ", plot_inside,  &
                                         vw_save(-grid % n_bnd_cells),         &
                                         f8, f9)
    end if

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
