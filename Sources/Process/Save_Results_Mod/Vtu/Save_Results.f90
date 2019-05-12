!==============================================================================!
  subroutine Save_Results(flow, turb, name_save)
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
                      kin_vis_t => r_cell_12
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  character(len=*)         :: name_save
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c, n, offset
  character(len=80)        :: name_out_8, name_out_9, store_name
!-----------------------------[Local parameters]-------------------------------!
  integer, parameter :: VTK_TETRA      = 10  ! cell shapes in VTK format
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

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  call Comm_Mod_Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  call Name_File(0, name_out_8, '.pvtu')
  call Name_File(this_proc, name_out_9, '.vtu')

  if(n_proc > 1 .and. this_proc .eq. 1) then
    open(8, file=name_out_8)
    print *, '# Creating file: ', trim(name_out_8)
  end if
  open(9, file=name_out_9)
  print *, '# Creating file: ', trim(name_out_9)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(8,'(a,a)') IN_0, '<?xml version="1.0"?>'
    write(8,'(a,a)') IN_0, '<VTKFile type="PUnstructuredGrid">'
    write(8,'(a,a)') IN_1, '<PUnstructuredGrid GhostLevel="0">'
  end if

  write(9,'(a,a)') IN_0, '<?xml version="1.0"?>'
  write(9,'(a,a)') IN_0, '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                         'byte_order="LittleEndian">'
  write(9,'(a,a)') IN_1, '<UnstructuredGrid>'

  write(9,'(a,a,i0.0,a,i0.0,a)')   &
               IN_2, '<Piece NumberOfPoints="', grid % n_nodes,      &
                          '" NumberOfCells ="', grid % n_cells -     &
                                                grid % comm % n_buff_cells, '">'

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  !-----------!
  !   Nodes   !
  !-----------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(8,'(a,a)') IN_3, '<PPoints>'
    write(8,'(a,a)') IN_4, '<PDataArray type="Float64" NumberOfComponents=' // &
                           '"3" format="ascii"/>'
    write(8,'(a,a)') IN_3, '</PPoints>'
  end if
  write(9,'(a,a)') IN_3, '<Points>'
  write(9,'(a,a)') IN_4, '<DataArray type="Float64" NumberOfComponents' //  &
                         '="3" format="ascii">'
  do n = 1, grid % n_nodes
    write(9, '(a,1pe16.6e4,1pe16.6e4,1pe16.6e4)')                &
               IN_5, grid % xn(n), grid % yn(n), grid % zn(n)
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Points>'

  !-----------!
  !   Cells   !
  !-----------!
  write(9,'(a,a)') IN_3, '<Cells>'

  ! First write all cells' nodes
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="connectivity"' //  &
                         ' format="ascii">'

  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % cells_n_nodes(c) .eq. 8) then
      write(9,'(a,8i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(4,c)-1, grid % cells_n(3,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1,   &
        grid % cells_n(8,c)-1, grid % cells_n(7,c)-1
    else if(grid % cells_n_nodes(c) .eq. 6) then
      write(9,'(a,6i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1,   &
        grid % cells_n(5,c)-1, grid % cells_n(6,c)-1
    else if(grid % cells_n_nodes(c) .eq. 4) then
      write(9,'(a,4i9)')                                &
        IN_5,                                           &
        grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
        grid % cells_n(3,c)-1, grid % cells_n(4,c)-1
    else if(grid % cells_n_nodes(c) .eq. 5) then
      write(9,'(a,5i9)')                                &
        IN_5,                                           &
        grid % cells_n(5,c)-1, grid % cells_n(1,c)-1,   &
        grid % cells_n(2,c)-1, grid % cells_n(4,c)-1,   &
        grid % cells_n(3,c)-1
    else
      print *, '# EERROR!  Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      call Comm_Mod_End
    end if
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' offsets
  write(9,'(a,a)') IN_4, '<DataArray type="Int64" Name="offsets" format="ascii">'
  offset = 0
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    offset = offset + grid % cells_n_nodes(c)
    write(9,'(a,i9)') IN_5, offset
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'

  ! Now write all cells' types
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="types" format="ascii">'
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % cells_n_nodes(c) .eq. 8) then
      write(9,'(a,i9)') IN_5, VTK_HEXAHEDRON
    else if(grid % cells_n_nodes(c) .eq. 6) then
      write(9,'(a,i9)') IN_5, VTK_WEDGE
    else if(grid % cells_n_nodes(c) .eq. 4) then
      write(9,'(a,i9)') IN_5, VTK_TETRA
    else if(grid % cells_n_nodes(c) .eq. 5) then
      write(9,'(a,i9)') IN_5, VTK_PYRAMID
    else
      print *, '# ERROR!  Unsupported cell type with ',  &
                  grid % cells_n_nodes(c), ' nodes.'
      print *, '# Exiting'
      call Comm_Mod_End
    end if
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'
  write(9,'(a,a)') IN_3, '</Cells>'

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(8,'(a,a)') IN_3, '<PCellData Scalars="scalars" vectors="velocity">'
  end if
  write(9,'(a,a)') IN_3, '<CellData Scalars="scalars" vectors="velocity">'

  !--------------------!
  !   Processor i.d.   !
  !--------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(8,'(a,a)') IN_3, '<PDataArray type="UInt8" Name="Processor"' //  &
                           ' format="ascii"/>'
  end if
  write(9,'(a,a)') IN_4, '<DataArray type="UInt8" Name="Processor"' //  &
                         ' format="ascii">'
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    write(9,'(a,i9)') IN_5, grid % comm % proces(c)
  end do
  write(9,'(a,a)') IN_4, '</DataArray>'

  !--------------!
  !   Velocity   !
  !--------------!
  call Save_Vector(grid, IN_4, IN_5, "Velocity",  &
                   flow % u % n(1), flow % v % n(1), flow % w % n(1))

  !--------------!
  !   Pressure   !
  !--------------!
  call Save_Scalar(grid, IN_4, IN_5, "Pressure", flow % p % n(1))

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(heat_transfer) then
    call Save_Scalar(grid, IN_4, IN_5, "Temperature", flow % t % n(1))
  end if

  !--------------------------!
  !   Turbulent quantities   !
  !--------------------------!

  ! Save kin and eps abd t2 for K_EPS (to improve)
  if(turbulence_model .eq. K_EPS                 .or.  &
     turbulence_model .eq. K_EPS_ZETA_F          .or.  &
     turbulence_model .eq. HYBRID_LES_RANS       .or.  &
     turbulence_model .eq. RSM_MANCEAU_HANJALIC  .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC  ) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergy",  &
                                           turb % kin % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentDissipation",    &
                                           turb % eps % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergyProduction", &
                                           p_kin(1))
      
    if (turbulence_model .eq. K_EPS .and. heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityT2",  &
                                             turb % t2 % n(1))
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentT2Production", &
                                             p_t2(1))
    end if

  end if

  ! Save zeta and f22
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells
      v2_calc(c) = turb % kin % n(c) * turb % zeta % n(c)
    end do
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityV2",   v2_calc (1))
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityZeta",  &
                                           turb % zeta % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22",   &
                                           turb % f22  % n(1))
    
    if (heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityT2",  &
                                             turb % t2 % n(1))
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentT2Production", &
                                             p_t2(1))
    end if
  
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22",   &
                                           turb % f22 % n(1))
  end if

  ! Save vis and vis_t
  if(turbulence_model .eq. DES_SPALART .or.  &
     turbulence_model .eq. SPALART_ALLMARAS) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentViscosity",  &
                                           turb % vis % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "VorticityMagnitude", flow % vort(1))
  end if
  if(turbulence_model .ne. NONE) then
    kin_vis_t(1:grid % n_cells) = vis_t(1:grid % n_cells)/viscosity
    call Save_Scalar(grid, IN_4, IN_5, "EddyOverMolecularViscosity", &
      kin_vis_t(1))
  end if

  ! Reynolds stress models
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXX", turb % uu % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYY", turb % vv % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressZZ", turb % ww % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXY", turb % uv % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXZ", turb % uw % n(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYZ", turb % vw % n(1))
  end if

  ! Statistics for large-scale simulations of turbulence
  if(turbulence_model .eq. LES_SMAGORINSKY    .or.  &
     turbulence_model .eq. LES_DYNAMIC        .or.  &
     turbulence_model .eq. LES_WALE           .or.  &
     turbulence_model .eq. DNS                .or.  &
     turbulence_model .eq. DES_SPALART        .or.  &
     turbulence_model .eq. HYBRID_LES_PRANDTL .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    call Save_Vector(grid, IN_4, IN_5, "MeanVelocity",  &
                                           turb % u_mean(1),  &
                                           turb % v_mean(1),  &
                                           turb % w_mean(1))
    do c = 1, grid % n_cells
      uu_save(c) = turb % uu_res(c) - turb % u_mean(c) * turb % u_mean(c)
      vv_save(c) = turb % vv_res(c) - turb % v_mean(c) * turb % v_mean(c)
      ww_save(c) = turb % ww_res(c) - turb % w_mean(c) * turb % w_mean(c)
      uv_save(c) = turb % uv_res(c) - turb % u_mean(c) * turb % v_mean(c)
      uw_save(c) = turb % uw_res(c) - turb % u_mean(c) * turb % w_mean(c)
      vw_save(c) = turb % vw_res(c) - turb % v_mean(c) * turb % w_mean(c)
    end do
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXX", uu_save(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYY", vv_save(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressZZ", ww_save(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXY", uv_save(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressXZ", uw_save(1))
    call Save_Scalar(grid, IN_4, IN_5, "ReynoldsStressYZ", vw_save(1))
    if(heat_transfer) then
      call Save_Scalar(grid, IN_4, IN_5, "TemperatureMean",  &
                                             turb % t_mean(1))
      do c = 1, grid % n_cells
        t2_save(c) = turb % t2_res(c) - turb % t_mean(c) * turb % t_mean(c)
        ut_save(c) = turb % ut_res(c) - turb % u_mean(c) * turb % t_mean(c)
        vt_save(c) = turb % vt_res(c) - turb % v_mean(c) * turb % t_mean(c)
        wt_save(c) = turb % wt_res(c) - turb % w_mean(c) * turb % t_mean(c)
      end do
      call Save_Scalar(grid, IN_4, IN_5, "TemperatureFluctuations",  &
                                 t2_save(1))
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxX", ut_save(1))
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxY", vt_save(1))
      call Save_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxZ", wt_save(1))
    end if
  end if

  ! Save y+ for all turbulence models
  if(turbulence_model .ne. NONE) then
    call Save_Scalar(grid, IN_4, IN_5, "TurbulentQuantityYplus", y_plus(1))
  end if

  ! Wall distance and delta, important for all models
  call Save_Scalar(grid, IN_4, IN_5, "WallDistance",  grid % wall_dist(1))
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaMax",  turb % h_max(1))
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaMin",  turb % h_min(1))
  call Save_Scalar(grid, IN_4, IN_5, "CellDeltaWall", turb % h_w  (1))

  !----------------------!
  !   Save user arrays   !
  !----------------------!
  ! call User_Mod_Save_Results(grid)

  !----------------------!
  !                      !
  !   End of cell data   !
  !                      !
  !----------------------!
  if(n_proc > 1 .and. this_proc .eq. 1) then
    write(8,'(a,a)') IN_3, '</PCellData>'
   end if
  write(9,'(a,a)') IN_3, '</CellData>'

  !------------!
  !            !
  !   Footer   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1) then
    do n = 1, n_proc
      call Name_File(n, name_out_9, '.vtu')
      write(8, '(a,a,a,a)') IN_2, '<Piece Source="', trim(name_out_9), '"/>'
    end do
    write(8, '(a,a)') IN_1, '</PUnstructuredGrid>'
    write(8, '(a,a)') IN_0, '</VTKFile>'
    close(8)
  end if
  write(9,'(a,a)') IN_2, '</Piece>'
  write(9,'(a,a)') IN_1, '</UnstructuredGrid>'
  write(9,'(a,a)') IN_0, '</VTKFile>'
  close(9)

  ! Restore the name
  problem_name = store_name

  call Cpu_Timer_Mod_Stop('Save_Vtu_Results')

  end subroutine
