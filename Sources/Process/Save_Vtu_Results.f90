!==============================================================================!
  subroutine Save_Results(flow, name_save)
!------------------------------------------------------------------------------!
!   Writes results in VTU file format (for VisIt and Paraview)                 !
!------------------------------------------------------------------------------!
!---------------------------------[Modules]------------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Field_Mod
  use Rans_Mod
  use Comm_Mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
  use User_Mod
  use Work_Mod, only: v2_calc   => r_cell_01,  &
                      uu_mean   => r_cell_02,  &
                      vv_mean   => r_cell_03,  &
                      ww_mean   => r_cell_04,  &
                      uv_mean   => r_cell_05,  &
                      uw_mean   => r_cell_06,  &
                      vw_mean   => r_cell_07,  &
                      tt_mean   => r_cell_08,  &
                      ut_mean   => r_cell_09,  &
                      vt_mean   => r_cell_10,  &
                      wt_mean   => r_cell_11,  &
                      kin_vis_t => r_cell_12
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type), target :: flow
  character(len=*)         :: name_save
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: phi
  integer                  :: c, n, offset, sc
  character(len=80)        :: name_out_8, name_out_9, store_name, name_mean
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
  call Save_Vtu_Vector(grid, IN_4, IN_5, "Velocity",  &
                       flow % u % n(1), flow % v % n(1), flow % w % n(1))

  !--------------!
  !   Pressure   !
  !--------------!
  call Save_Vtu_Scalar(grid, IN_4, IN_5, "Pressure", flow % p % n(1))

  !-----------------!
  !   Temperature   !
  !-----------------!
  if(heat_transfer) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "Temperature", flow % t % n(1))
  end if

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, flow % n_scalars
    phi => flow % scalar(sc)
    call Save_Vtu_Scalar(grid, IN_4, IN_5, phi % name, phi % n(1))
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
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergy", kin % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentDissipation",   eps % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentKineticEnergyProduction", &
                                           p_kin(1))
  end if

  ! Save zeta and f22
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells
      v2_calc(c) = kin % n(c) * zeta % n(c)
    end do
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityV2",   v2_calc (1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityZeta", zeta % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22",  f22  % n(1))
  end if
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityF22",  f22  % n(1))
  end if

  if( (turbulence_model .eq. K_EPS_ZETA_F .and. heat_transfer) .or. &
      (turbulence_model .eq. HYBRID_LES_RANS .and. heat_transfer) ) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityT2", t2 % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentT2Production", p_t2(1))
    if(buoyancy) &
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentBuoyancyProduction", g_buoy(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxX", ut % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxY", vt % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxZ", wt % n(1))
  end if
  ! Save vis and vis_t
  if(turbulence_model .eq. DES_SPALART .or.  &
     turbulence_model .eq. SPALART_ALLMARAS) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentViscosity", vis % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "VorticityMagnitude", vort(1))
  end if
  if(turbulence_model .ne. NONE) then                      
    kin_vis_t(1:grid % n_cells) = vis_t(1:grid % n_cells)/viscosity
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "EddyOverMolecularViscosity", &
      kin_vis_t(1))
  end if

  ! Reynolds stress models
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXX", uu % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressYY", vv % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressZZ", ww % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXY", uv % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXZ", uw % n(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressYZ", vw % n(1))
  end if

  ! Statistics for large-scale simulations of turbulence
  if(turbulence_statistics) then
    call Save_Vtu_Vector(grid, IN_4, IN_5, "MeanVelocity",  &
                                           flow % u % mean(1),  &
                                           flow % v % mean(1),  &
                                           flow % w % mean(1))
    do c = 1, grid % n_cells
      uu_mean(c) = uu % mean(c) - flow % u % mean(c) * flow % u % mean(c)
      vv_mean(c) = vv % mean(c) - flow % v % mean(c) * flow % v % mean(c)
      ww_mean(c) = ww % mean(c) - flow % w % mean(c) * flow % w % mean(c)
      uv_mean(c) = uv % mean(c) - flow % u % mean(c) * flow % v % mean(c)
      uw_mean(c) = uw % mean(c) - flow % u % mean(c) * flow % w % mean(c)
      vw_mean(c) = vw % mean(c) - flow % v % mean(c) * flow % w % mean(c)
    end do
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXX", uu_mean(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressYY", vv_mean(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressZZ", ww_mean(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXY", uv_mean(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressXZ", uw_mean(1))
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "ReynoldsStressYZ", vw_mean(1))

    if(heat_transfer) then
      call Save_Vtu_Scalar(grid, IN_4, IN_5, "TemperatureMean",  &
                                             flow % t % mean(1))
      do c = 1, grid % n_cells
        tt_mean(c) = tt % mean(c) - flow % t % mean(c) * flow % t % mean(c)
        ut_mean(c) = ut % mean(c) - flow % u % mean(c) * flow % t % mean(c)
        vt_mean(c) = vt % mean(c) - flow % v % mean(c) * flow % t % mean(c)
        wt_mean(c) = wt % mean(c) - flow % w % mean(c) * flow % t % mean(c)
      end do
      call Save_Vtu_Scalar(grid, IN_4, IN_5, "TemperatureFluctuations",  &
                                 tt_mean(1))
      call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxX", ut_mean(1))
      call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxY", vt_mean(1))
      call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentHeatFluxZ", wt_mean(1))
    end if

    if(flow % n_scalars > 0) then
      do sc = 1, flow % n_scalars
        phi => flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        call Save_Vtu_Scalar(grid, IN_4, IN_5, name_mean, phi % mean(1))
      end do
    end if
  end if

  ! Save y+ for all turbulence models
  if(turbulence_model .ne. NONE) then
    call Save_Vtu_Scalar(grid, IN_4, IN_5, "TurbulentQuantityYplus", y_plus(1))
  end if

  ! Wall distance and delta, important for all models
  call Save_Vtu_Scalar(grid, IN_4, IN_5, "WallDistance", grid % wall_dist(1))
  call Save_Vtu_Scalar(grid, IN_4, IN_5, "CellDelta",    grid % delta(1))

  !----------------------!
  !   Save user arrays   !
  !----------------------!
  call User_Mod_Save_Vtu_Results(grid)

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

  end subroutine
