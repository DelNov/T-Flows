!==============================================================================!
  subroutine Save_Results(flow, turb, mult, swarm, ts, plot_inside)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
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
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: ts           ! time step
  logical                       :: plot_inside  ! plot results inside?
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: phi
  integer                  :: c, n, s, offset, sc, f8, f9, ua
  character(len=80)        :: name_out_8, name_out_9, name_mean, a_name
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: VTK_TRIANGLE   =  5  ! cell shapes in VTK format
  integer, parameter :: VTK_QUAD       =  9
  integer, parameter :: VTK_TETRA      = 10
  integer, parameter :: VTK_HEXAHEDRON = 12
  integer, parameter :: VTK_WEDGE      = 13
  integer, parameter :: VTK_PYRAMID    = 14
  character(len= 0)  :: IN_0 = ''           ! indentation levels
  character(len= 2)  :: IN_1 = '  '
  character(len= 4)  :: IN_2 = '    '
  character(len= 6)  :: IN_3 = '      '
  character(len= 8)  :: IN_4 = '        '
  character(len=10)  :: IN_5 = '          '
!==============================================================================!

  call Cpu_Timer_Mod_Start('Save_Vtu_Results')

  ! Take aliases
  grid => flow % pnt_grid

  call Comm_Mod_Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File_Mod_Set_Name(name_out_8,             &
                           time_step=ts,           &
                           extension='.pvtu')
    call File_Mod_Set_Name(name_out_9,             &
                           processor=this_proc,    &
                           time_step=ts,           &
                           extension='.vtu')
  else
    call File_Mod_Set_Name(name_out_8,             &
                           time_step=ts,           &
                           appendix ='-bnd',       &
                           extension='.pvtu')
    call File_Mod_Set_Name(name_out_9,             &
                           processor=this_proc,    &
                           time_step=ts,           &
                           appendix ='-bnd',       &
                           extension='.vtu')
  end if

  if(n_proc > 1 .and. this_proc .eq. 1) then
    call File_Mod_Open_File_For_Writing(name_out_8, f8)
  end if
  call File_Mod_Open_File_For_Writing(name_out_9, f9)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(f8,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(f8,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'
  end if

  write(f9,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(f9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' // &
                          'byte_order="LittleEndian">'
  write(f9,'(a,a)') IN_1, '<UnstructuredGrid>'

  if(plot_inside) then
    write(f9,'(a,a,i0.0,a,i0.0,a)')   &
                 IN_2, '<Piece NumberOfPoints="', grid % n_nodes,              &
                            '" NumberOfCells ="', grid % n_cells               &
                                                - grid % comm % n_buff_cells,  &
                                                '">'
  else
    write(f9,'(a,a,i0.0,a,i0.0,a)')   &
                 IN_2, '<Piece NumberOfPoints="', grid % n_nodes,              &
                            '" NumberOfCells ="', grid % n_bnd_cells, '">'
  end if

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  !-----------!
  !   Nodes   !
  !-----------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8,'(a,a)') IN_3, '<PPoints>'
    write(f8,'(a,a)') IN_4, '<PDataArray type="Float64" NumberOfComponents='//&
                            '"3" format="ascii"/>'
    write(f8,'(a,a)') IN_3, '</PPoints>'
  end if
  write(f9,'(a,a)') IN_3, '<Points>'
  write(f9,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                          '="3" format="ascii">'
  do n = 1, grid % n_nodes
    write(f9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
                IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(f9,'(a,a)') IN_4, '</DataArray>'
  write(f9,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Cells   !
  !-----------!
  write(f9,'(a,a)') IN_3, '<Cells>'

  ! First write all cells' nodes
  write(f9,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                          ' format="ascii">'

  if(plot_inside) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      if(grid % cells_n_nodes(c) .eq. 8) then
        write(f9,'(a,8i9)')                                &
           IN_5,                                           &
           grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
           grid % cells_n(4,c)-1, grid % cells_n(3,c)-1,   &
           grid % cells_n(5,c)-1, grid % cells_n(6,c)-1,   &
           grid % cells_n(8,c)-1, grid % cells_n(7,c)-1
      else if(grid % cells_n_nodes(c) .eq. 6) then
        write(f9,'(a,6i9)')                                &
           IN_5,                                           &
           grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
           grid % cells_n(3,c)-1, grid % cells_n(4,c)-1,   &
           grid % cells_n(5,c)-1, grid % cells_n(6,c)-1
      else if(grid % cells_n_nodes(c) .eq. 4) then
        write(f9,'(a,4i9)')                                &
           IN_5,                                           &
           grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
           grid % cells_n(3,c)-1, grid % cells_n(4,c)-1
      else if(grid % cells_n_nodes(c) .eq. 5) then
        write(f9,'(a,5i9)')                                &
           IN_5,                                           &
           grid % cells_n(5,c)-1, grid % cells_n(1,c)-1,   &
           grid % cells_n(2,c)-1, grid % cells_n(4,c)-1,   &
           grid % cells_n(3,c)-1
      else
        print *, '# ERROR!  Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        call Comm_Mod_End
      end if
    end do
  else  ! plot only boundary
    do s = 1, grid % n_faces
      if( grid % faces_c(2,s) < 0 ) then
        if(grid % faces_n_nodes(s) .eq. 4) then
          write(f9,'(a,4I9)')                               &
             IN_5,                                          &
             grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
             grid % faces_n(3,s)-1, grid % faces_n(4,s)-1
        else if(grid % faces_n_nodes(s) .eq. 3) then
          write(f9,'(a,3I9)')                               &
             IN_5,                                          &
             grid % faces_n(1,s)-1, grid % faces_n(2,s)-1,  &
             grid % faces_n(3,s)-1
        else
          print *, '# ERROR!  Unsupported face type ',      &
                    grid % faces_n_nodes(s), ' nodes.'
          print *, '# Exiting'
          call Comm_Mod_End
        end if
      end if
    end do
  end if

  write(f9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(f9,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets"' //  &
                          ' format="ascii">'
  offset = 0

  if(plot_inside) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      offset = offset + grid % cells_n_nodes(c)
      write(f9,'(a,i9)') IN_5, offset
    end do
  else  ! plot only boundary
    do s = 1, grid % n_faces
      if( grid % faces_c(2,s) < 0 ) then
        offset = offset + grid % faces_n_nodes(s)
        write(f9,'(a,i9)') IN_5, offset
      end if
    end do
  end if
  write(f9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(f9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="types" format="ascii">'

  if(plot_inside) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      if(grid % cells_n_nodes(c) .eq. 8) then
        write(f9,'(a,i9)') IN_5, VTK_HEXAHEDRON
      else if(grid % cells_n_nodes(c) .eq. 6) then
        write(f9,'(a,i9)') IN_5, VTK_WEDGE
      else if(grid % cells_n_nodes(c) .eq. 4) then
        write(f9,'(a,i9)') IN_5, VTK_TETRA
      else if(grid % cells_n_nodes(c) .eq. 5) then
        write(f9,'(a,i9)') IN_5, VTK_PYRAMID
      else
        print *, '# ERROR!  Unsupported cell type with ',  &
                    grid % cells_n_nodes(c), ' nodes.'
        print *, '# Exiting'
        call Comm_Mod_End
      end if
    end do
  else  ! plot only boundary
    do s = 1, grid % n_faces
      if( grid % faces_c(2,s) < 0 ) then
        if(grid % faces_n_nodes(s) .eq. 4) write(f9,'(a,i9)') IN_5, VTK_QUAD
        if(grid % faces_n_nodes(s) .eq. 3) write(f9,'(a,i9)') IN_5, VTK_TRIANGLE
      end if
    end do
  end if

  write(f9,'(a,a)') IN_4, '</DataArray>'
  write(f9,'(a,a)') IN_3, '</Cells>'

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8,'(a,a)') IN_3, '<PCellData Scalars="scalars" vectors="velocity">'
  end if
  write(f9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  !--------------------!
  !   Processor i.d.   !
  !--------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8,'(a,a)') IN_3, '<PDataArray type="UInt8" Name="Processor"' //  &
                            ' format="ascii"/>'
  end if
  write(f9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="Processor"' //  &
                          ' format="ascii">'
  if(plot_inside) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      write(f9,'(a,i9)') IN_5, grid % comm % cell_proc(c)
    end do
  else  ! plot only boundary
    do s = 1, grid % n_faces
      if( grid % faces_c(2,s) < 0 ) then
        write(f9,'(a,i9)') IN_5, grid % comm % cell_proc( grid % faces_c(1,s) )
      end if
    end do
  end if
  write(f9,'(a,a)') IN_4, '</DataArray>'

  !--------------!
  !   Velocity   !
  !--------------!
  call Save_Vector(grid, IN_4, IN_5, "Velocity", plot_inside,            &
                                     flow % u % n(-grid % n_bnd_cells),  &
                                     flow % v % n(-grid % n_bnd_cells),  &
                                     flow % w % n(-grid % n_bnd_cells),  &
                                     f8, f9)

  !--------------!
  !   Pressure   !
  !--------------!
  call Save_Scalar(grid, IN_4, IN_5, "Pressure", plot_inside,            &
                                     flow % p % n(-grid % n_bnd_cells),  &
                                     f8, f9)

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(heat_transfer) then
    call Save_Scalar(grid, IN_4, IN_5, "Temperature", plot_inside,         &
                                       flow % t % n(-grid % n_bnd_cells),  &
                                       f8, f9)
  end if

  !---------------------!
  !   Volume Fraction   !
  !---------------------!
  if(multiphase_model .eq. VOLUME_OF_FLUID) then
    call Save_Scalar(grid, IN_4, IN_5, "VolumeFraction", plot_inside,        &
                                       mult % vof % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    if (mult % d_func) then
      call Save_Scalar(grid, IN_4, IN_5, "DistanceFunction", plot_inside,    &
                                         mult % dist_func                    &
                                         % n(-grid % n_bnd_cells),           &
                                         f8, f9)
    end if
  end if

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, flow % n_scalars
    phi => flow % scalar(sc)
    call Save_Scalar(grid, IN_4, IN_5, phi % name, plot_inside,  &
                                       phi % n(-grid % n_bnd_cells), f8, f9)
  end do

  !--------------------------!
  !   Turbulent quantities   !
  !--------------------------!

  ! Save kin and eps
  if(turbulence_model .eq. K_EPS                 .or.  &
     turbulence_model .eq. K_EPS_ZETA_F          .or.  &
     turbulence_model .eq. HYBRID_LES_RANS       .or.  &
     turbulence_model .eq. RSM_MANCEAU_HANJALIC  .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC  ) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergy", plot_inside,  &
                                       turb % kin % n(-grid % n_bnd_cells),    &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentDissipation", plot_inside,    &
                                       turb % eps % n(-grid % n_bnd_cells),    &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergyProduction",     &
                                       plot_inside,                            &
                                       turb % p_kin(-grid % n_bnd_cells),      &
                                       f8, f9)
  end if

  ! Save zeta and f22
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    v2_calc(:) = 0.0
    do c = 1, grid % n_cells
      v2_calc(c) = turb % kin % n(c) * turb % zeta % n(c)
    end do
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityV2", plot_inside,     &
                                       v2_calc (-grid % n_bnd_cells),          &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityZeta", plot_inside,   &
                                       turb % zeta % n(-grid % n_bnd_cells),   &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22", plot_inside,    &
                                       turb % f22  % n(-grid % n_bnd_cells),   &
                                       f8, f9)
    if (heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityT2", plot_inside,   &
                                         turb % t2 % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentT2Production", plot_inside, &
                                         turb % p_t2(-grid % n_bnd_cells),     &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxX", plot_inside,    &
                                         turb % ut % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxY", plot_inside,    &
                                         turb % vt % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxZ", plot_inside,    &
                                         turb % wt % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulenQuantityAlphaL",             &
                                         plot_inside,                          &
                                         turb % alpha_l(-grid % n_bnd_cells),  &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulenQuantityAlphaU",             &
                                         plot_inside,                          &
                                         turb % alpha_u(-grid % n_bnd_cells),  &
                                         f8, f9)
    end if
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22", plot_inside,  &
                                       turb % f22 % n(-grid % n_bnd_cells),  &
                                       f8, f9)
  end if

  ! Save vis and vis_t
  if(turbulence_model .eq. DES_SPALART .or.  &
     turbulence_model .eq. SPALART_ALLMARAS) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentViscosity", plot_inside,    &
                                       turb % vis % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "VorticityMagnitude", plot_inside,    &
                                       flow % vort(-grid % n_bnd_cells),     &
                                       f8, f9)
  end if
  kin_vis_t(:) = 0.0
  if(turbulence_model .ne. NO_TURBULENCE) then
    kin_vis_t   (-grid % n_bnd_cells:grid % n_cells) =  &
    turb % vis_t(-grid % n_bnd_cells:grid % n_cells) /  &
       flow % viscosity(-grid % n_bnd_cells:grid % n_cells)
    call Save_Scalar(grid, IN_4, IN_5, "EddyOverMolecularViscosity",        &
                                       plot_inside,                         &
                                       kin_vis_t(-grid % n_bnd_cells), f8, f9)
  end if

  ! Reynolds stress models
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXX", plot_inside,     &
                                       turb % uu % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYY", plot_inside,     &
                                       turb % vv % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressZZ", plot_inside,     &
                                       turb % ww % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXY", plot_inside,     &
                                       turb % uv % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXZ", plot_inside,     &
                                       turb % uw % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYZ", plot_inside,     &
                                       turb % vw % n(-grid % n_bnd_cells),  &
                                       f8, f9)
    if(heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxX", plot_inside,    &
                                         turb % ut % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxY", plot_inside,    &
                                         turb % vt % n(-grid % n_bnd_cells),   &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxZ", plot_inside,    &
                                         turb % wt % n(-grid % n_bnd_cells),   &
                                         f8, f9)
    end if
  end if

  ! Statistics for large-scale simulations of turbulence
  if(turbulence_statistics) then
    call Save_Vector(grid, IN_4, IN_5, "MeanVelocity", plot_inside,          &
                                       turb % u_mean(-grid % n_bnd_cells),   &
                                       turb % v_mean(-grid % n_bnd_cells),   &
                                       turb % w_mean(-grid % n_bnd_cells),   &
                                       f8, f9)
    uu_save(:) = 0.0
    vv_save(:) = 0.0
    ww_save(:) = 0.0
    uv_save(:) = 0.0
    uw_save(:) = 0.0
    vw_save(:) = 0.0
    do c = 1, grid % n_cells
      uu_save(c) = turb % uu_res(c) - turb % u_mean(c) * turb % u_mean(c)
      vv_save(c) = turb % vv_res(c) - turb % v_mean(c) * turb % v_mean(c)
      ww_save(c) = turb % ww_res(c) - turb % w_mean(c) * turb % w_mean(c)
      uv_save(c) = turb % uv_res(c) - turb % u_mean(c) * turb % v_mean(c)
      uw_save(c) = turb % uw_res(c) - turb % u_mean(c) * turb % w_mean(c)
      vw_save(c) = turb % vw_res(c) - turb % v_mean(c) * turb % w_mean(c)
    end do
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressXX", plot_inside,  &
                                       uu_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressYY", plot_inside,  &
                                       vv_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressZZ", plot_inside,  &
                                       ww_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressXY", plot_inside,  &
                                       uv_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressXZ", plot_inside,  &
                                       uw_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    call Save_Scalar(grid, IN_4, IN_5, "MeanReynoldsStressYZ", plot_inside,  &
                                       vw_save(-grid % n_bnd_cells),         &
                                       f8, f9)
    if(heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "MeanTemperature", plot_inside,      &
                                         turb % t_mean(-grid % n_bnd_cells),  &
                                         f8, f9)
      t2_save(:) = 0.0
      ut_save(:) = 0.0
      vt_save(:) = 0.0
      wt_save(:) = 0.0
      do c = 1, grid % n_cells
        t2_save(c) = turb % t2_res(c) - turb % t_mean(c) * turb % t_mean(c)
        ut_save(c) = turb % ut_res(c) - turb % u_mean(c) * turb % t_mean(c)
        vt_save(c) = turb % vt_res(c) - turb % v_mean(c) * turb % t_mean(c)
        wt_save(c) = turb % wt_res(c) - turb % w_mean(c) * turb % t_mean(c)
      end do
      call Save_Scalar(grid, IN_4, IN_5, "MeanTurbulentQuantityT2",     &
                                         plot_inside,                   &
                                         t2_save(-grid % n_bnd_cells),  &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanTurbulentHeatFluxX",      &
                                         plot_inside,                   &
                                         ut_save(-grid % n_bnd_cells),  &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanTurbulentHeatFluxY",      &
                                         plot_inside,                   &
                                         vt_save(-grid % n_bnd_cells),  &
                                         f8, f9)
      call Save_Scalar(grid, IN_4, IN_5, "MeanTurbulentHeatFluxZ",      &
                                         plot_inside,                   &
                                         wt_save(-grid % n_bnd_cells),  &
                                         f8, f9)
    end if

    ! Scalars
    do sc = 1, flow % n_scalars
      phi => flow % scalar(sc)
      name_mean = 'Mean'
      name_mean(5:8) = phi % name
      do c = 1, grid % n_cells
        phi_save(c) = turb % scalar_mean(sc, c)
      end do
      call Save_Scalar(grid, IN_4, IN_5, name_mean, plot_inside,  &
                       phi_save(-grid % n_bnd_cells), f8, f9)
    end do
  end if

  ! Save y+ for all turbulence models
  if(turbulence_model .ne. NO_TURBULENCE) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityYplus",            &
                                       plot_inside,                         &
                                       turb % y_plus(-grid % n_bnd_cells),  &
                                       f8, f9)
  end if

  ! Wall distance and delta, important for all models
  call Save_Scalar(grid, IN_4, IN_5, "WallDistance", plot_inside,            &
                                     grid % wall_dist(-grid % n_bnd_cells),  &
                                     f8, f9)
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaMax", plot_inside,        &
                                     turb % h_max(-grid % n_bnd_cells),  &
                                     f8, f9)
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaMin", plot_inside,        &
                                     turb % h_min(-grid % n_bnd_cells),  &
                                     f8, f9)
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaWall", plot_inside,       &
                                     turb % h_w  (-grid % n_bnd_cells),  &
                                     f8, f9)

  !----------------------!
  !   Save user arrays   !
  !----------------------!
  do ua = 1, n_user_arrays

    a_name = 'A_00'
    write(a_name(3:4), '(I2.2)') ua
    call Save_Scalar(grid, IN_4, IN_5, a_name,                              &
                                       plot_inside,                         &
                                       user_array(ua,-grid % n_bnd_cells),  &
                                       f8, f9)
  end do

  !----------------------!
  !                      !
  !   End of cell data   !
  !                      !
  !----------------------!
  if(n_proc > 1 .and. this_proc .eq. 1) then
    write(f8,'(a,a)') IN_3, '</PCellData>'
   end if
  write(f9,'(a,a)') IN_3, '</CellData>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1) then
    do n = 1, n_proc
      if(plot_inside) then
        call File_Mod_Set_Name(name_out_9,        &
                               processor=n,       &
                               time_step=ts,      &
                               extension='.vtu')
      else
        call File_Mod_Set_Name(name_out_9,        &
                               processor=n,       &
                               time_step=ts,      &
                               appendix ='-bnd',  &
                               extension='.vtu')
      end if
      write(f8, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out_9), '"/>'
    end do
    write(f8, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(f8, '(a,a)') IN_0, '</VTKFile>'
    close(f8)
  end if
  write(f9,'(a,a)') IN_2, '</Piece>'
  write(f9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(f9,'(a,a)') IN_0, '</VTKFile>'
  close(f9)

  call Cpu_Timer_Mod_Stop('Save_Vtu_Results')

  end subroutine
