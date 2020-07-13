!==============================================================================!
  subroutine Save_Results(flow, turb, mult, swarm, ts, plot_inside, domain)
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
                      phi_save  => r_cell_13,  &
                      q_save    => r_cell_14,  &
                      int_save  => i_cell_01,  &
                      type_save => i_cell_02,  &  ! cell type save array
                      offs_save => i_cell_03      ! cell offsets save array
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: ts           ! time step
  logical                       :: plot_inside  ! plot results inside?
  integer,             optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: phi
  integer(4)               :: data_size
  integer                  :: data_offset, cell_offset
  integer                  :: c, n, n_conns, sc, f8, f9, ua, run, c2
  character(len=80)        :: name_out_8, name_out_9, name_mean, a_name
  character(len=80)        :: str1, str2
  integer, parameter       :: IP=8, RP=8, SP=4
!==============================================================================!

  call Cpu_Timer_Mod_Start('Save_Vtu_Results')

  ! Take aliases
  grid => flow % pnt_grid

  ! Count connections in this grid, you will need it later
  n_conns = 0
  if(plot_inside) then
    do c = 1, grid % n_cells
      n_conns = n_conns + grid % cells_n_nodes(c)
    end do
  else
    do c2 = -grid % n_bnd_cells, -1
      n_conns = n_conns + grid % cells_n_nodes(c2)
    end do
  end if

  call Comm_Mod_Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File_Mod_Set_Name(name_out_8,             &
                           time_step=ts,           &
                           extension='.pvtu',      &
                           domain=domain)
    call File_Mod_Set_Name(name_out_9,             &
                           processor=this_proc,    &
                           time_step=ts,           &
                           extension='.vtu',       &
                           domain=domain)
  else
    call File_Mod_Set_Name(name_out_8,             &
                           time_step=ts,           &
                           appendix ='-bnd',       &
                           extension='.pvtu',      &
                           domain=domain)
    call File_Mod_Set_Name(name_out_9,             &
                           processor=this_proc,    &
                           time_step=ts,           &
                           appendix ='-bnd',       &
                           extension='.vtu',       &
                           domain=domain)
  end if

  if(n_proc > 1 .and. this_proc .eq. 1) then
    call File_Mod_Open_File_For_Writing_Binary(name_out_8, f8)
  end if
  call File_Mod_Open_File_For_Writing_Binary(name_out_9, f9)

  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_0 // '<?xml version="1.0"?>'              // LF
    write(f8) IN_0 // '<VTKFile type="PUnstructuredGrid">' // LF
    write(f8) IN_1 // '<PUnstructuredGrid GhostLevel="1">' // LF
  end if

  write(f9) IN_0 // '<?xml version="1.0"?>'                           // LF
  write(f9) IN_0 // '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                    'byte_order="LittleEndian">'                      // LF
  write(f9) IN_1 // '<UnstructuredGrid>'                              // LF

  write(str1,'(i0.0)') grid % n_nodes
  if(plot_inside) then
    write(str2,'(i0.0)') grid % n_cells
  else
    if(grid % n_bnd_cells .eq. 0) then
      write(str2,'(i1)')   grid % n_bnd_cells  ! i0.0 doesn't work for zero :-/
    else
      write(str2,'(i0.0)') grid % n_bnd_cells
    end if
  end if
  write(f9) IN_2 // '<Piece NumberOfPoints="' // trim(str1)    //  &
                    '" NumberOfCells ="' // trim(str2) // '">' // LF

  !----------!
  !          !
  !   Grid   !
  !          !
  !----------!

  data_offset = 0

  !-----------!
  !   Nodes   !
  !-----------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_3 // '<PPoints>' // LF
    write(f8) IN_4 // '<PDataArray type="Float64"' //  &
                      ' NumberOfComponents="3"/>'  // LF
    write(f8) IN_3 // '</PPoints>' // LF
  end if

  write(str1, '(i1)') data_offset
  write(f9) IN_3 // '<Points>'                       // LF
  write(f9) IN_4 // '<DataArray type="Float64"'      //  &
                    ' NumberOfComponents="3"'        //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(f9) IN_4 // '</DataArray>' // LF
  write(f9) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + grid % n_nodes * RP * 3

  !-----------!
  !   Cells   !
  !-----------!
  write(f9) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(f9) IN_4 // '<DataArray type="Int64"'        //  &
                    ' Name="connectivity"'           //  &
                    ' format="appended"'             //  &
                    ' offset="' // trim(str1) //'">' // LF
  write(f9) IN_4 // '</DataArray>' // LF
  data_offset = data_offset + SP + n_conns * IP  ! prepare for next

  ! Fill up an array with cell offsets and save the header only
  cell_offset = 0
  do c = 1, grid % n_cells
    cell_offset = cell_offset + grid % cells_n_nodes(c)
    offs_save(c) = cell_offset
  end do
  cell_offset = 0
  do c2 = -grid % n_bnd_cells, -1
    cell_offset = cell_offset + grid % cells_n_nodes(c2)
    offs_save(c2) = cell_offset
  end do
  call Save_Scalar_Int(grid, "offsets", plot_inside,           &
                              offs_save(-grid % n_bnd_cells),  &
                              f8, f9, data_offset, 1)  ! 1 => header only

  ! Fill up an array with cell types and save the header only
  do c = 1, grid % n_cells
    if(grid % cells_n_nodes(c) .eq. 8) type_save(c) = VTK_HEXAHEDRON
    if(grid % cells_n_nodes(c) .eq. 6) type_save(c) = VTK_WEDGE
    if(grid % cells_n_nodes(c) .eq. 4) type_save(c) = VTK_TETRA
    if(grid % cells_n_nodes(c) .eq. 5) type_save(c) = VTK_PYRAMID
  end do
  do c2 = -grid % n_bnd_cells, -1
    if(grid % cells_n_nodes(c2) .eq. 4) type_save(c2) = VTK_QUAD
    if(grid % cells_n_nodes(c2) .eq. 3) type_save(c2) = VTK_TRIANGLE
  end do
  call Save_Scalar_Int(grid, "types", plot_inside,             &
                              type_save(-grid % n_bnd_cells),  &
                              f8, f9, data_offset, 1)  ! 1 => header only

  write(f9) IN_3 // '</Cells>'     // LF

  !---------------------------------!
  !                                 !
  !   Results and other cell data   !
  !                                 !
  !---------------------------------!
  if(n_proc > 1 .and. this_proc .eq. 1)  then
    write(f8) IN_3 // '<PCellData Scalars="scalars" vectors="velocity">' // LF
  end if
  write(f9) IN_3 // '<CellData Scalars="scalars" vectors="velocity">' // LF

  !----------------!
  !                !
  !   Two sweeps   !
  !                !
  !----------------!
  do run = 1, 2

    !------------------------------------------!
    !   Save remnants of the grid definition   !
    !------------------------------------------!
    if(run .eq. 2) then

      ! Save the nodes' coordinates
      data_size = grid % n_nodes * RP * 3
      write(f9) data_size
      do n = 1, grid % n_nodes
        write(f9) grid % xn(n), grid % yn(n), grid % zn(n)
      end do

      ! Save connections
      data_size = n_conns * IP
      write(f9) data_size
      if(plot_inside) then
        do c = 1, grid % n_cells
          if(grid % cells_n_nodes(c) .eq. 8) then
            write(f9) grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
                      grid % cells_n(4,c)-1, grid % cells_n(3,c)-1,   &
                      grid % cells_n(5,c)-1, grid % cells_n(6,c)-1,   &
                      grid % cells_n(8,c)-1, grid % cells_n(7,c)-1
          else if(grid % cells_n_nodes(c) .eq. 6) then
            write(f9) grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
                      grid % cells_n(3,c)-1, grid % cells_n(4,c)-1,   &
                      grid % cells_n(5,c)-1, grid % cells_n(6,c)-1
          else if(grid % cells_n_nodes(c) .eq. 4) then
            write(f9) grid % cells_n(1,c)-1, grid % cells_n(2,c)-1,   &
                      grid % cells_n(3,c)-1, grid % cells_n(4,c)-1
          else if(grid % cells_n_nodes(c) .eq. 5) then
            write(f9) grid % cells_n(5,c)-1, grid % cells_n(1,c)-1,   &
                      grid % cells_n(2,c)-1, grid % cells_n(4,c)-1,   &
                      grid % cells_n(3,c)-1
          end if
        end do
      else  ! plot only boundary
        do c2 = -grid % n_bnd_cells, -1
          if(grid % cells_n_nodes(c2) .eq. 4) then
            write(f9) grid % cells_n(1,c2)-1, grid % cells_n(2,c2)-1,  &
                      grid % cells_n(3,c2)-1, grid % cells_n(4,c2)-1
          else if(grid % cells_n_nodes(c2) .eq. 3) then
            write(f9) grid % cells_n(1,c2)-1, grid % cells_n(2,c2)-1,  &
                      grid % cells_n(3,c2)-1
          end if
        end do
      end if

      ! Save cell offsets
      call Save_Scalar_Int(grid, "offsets", plot_inside,           &
                                  offs_save(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      ! Save cell types
      call Save_Scalar_Int(grid, "types", plot_inside,             &
                                  type_save(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    !--------------------!
    !   Processor i.d.   !
    !--------------------!
    do c = 1, grid % n_cells
      int_save(c) = grid % comm % cell_proc(c)
    end do
    do c2 = -grid % n_bnd_cells, -1
      int_save(c2) = grid % comm % cell_proc(c2)
    end do
    call Save_Scalar_Int(grid, "Processor", plot_inside,        &
                                int_save(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(-grid % n_bnd_cells:grid % n_cells) = domain
      call Save_Scalar_Int(grid, "Domain", plot_inside,           &
                                  int_save(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    !--------------!
    !   Velocity   !
    !--------------!
    call Save_Vector_Real(grid, "Velocity", plot_inside,            &
                                flow % u % n(-grid % n_bnd_cells),  &
                                flow % v % n(-grid % n_bnd_cells),  &
                                flow % w % n(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !---------------!
    !   Potential   !
    !---------------!
    call Save_Scalar_Real(grid, "Potential", plot_inside,             &
                                flow % pot % n(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !--------------!
    !   Pressure   !
    !--------------!
    call Save_Scalar_Real(grid, "Pressure", plot_inside,            &
                                flow % p % n(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !-----------------!
    !   Temperature   !
    !-----------------!
    if(heat_transfer) then
      call Save_Scalar_Real(grid, "Temperature", plot_inside,         &
                                  flow % t % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    !---------------------!
    !   Volume fraction   !
    !---------------------!
    if(mult % model .eq. VOLUME_OF_FLUID) then
      call Save_Scalar_Real(grid, "VolumeFraction", plot_inside,        &
                                  mult % vof % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "Curvature ", plot_inside,        &
                                  mult % vof % oo(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      if (allocated(mult % flux_rate)) then
        call Save_Scalar_Real(grid, "Flux_rate ", plot_inside,        &
                                    mult % flux_rate(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
      end if

      if (mult % d_func) then
        call Save_Scalar_Real(grid, "DistanceFunction", plot_inside,    &
                                    mult % dist_func                    &
                                    % n(-grid % n_bnd_cells),           &
                                    f8, f9, data_offset, run)
      end if
    end if

    !---------------------------------------!
    !   Number of impacts and reflections   !
    !---------------------------------------!
    if(mult % model .eq. LAGRANGIAN_PARTICLES .and. .not. plot_inside) then
      call Save_Scalar_Real(grid, "ParticlesReflected", plot_inside,         &
                                  swarm % n_reflected(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ParticlesDeposited", plot_inside,         &
                                  swarm % n_deposited(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, flow % n_scalars
      phi => flow % scalar(sc)
      call Save_Scalar_Real(grid, phi % name, plot_inside,       &
                                  phi % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end do

    !-----------------!
    !   Q-criterion   !
    !-----------------!
    q_save(:) = 0.0
    do c = 1, grid % n_cells
      q_save(c) = (flow % vort(c)**2 - flow % shear(c)**2)/4.
    end do
    call Save_Scalar_Real(grid, "QCriterion", plot_inside,   &
                                q_save(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !--------------------------!
    !   Turbulent quantities   !
    !--------------------------!

    ! Save kin and eps
    if(turb % model .eq. K_EPS                 .or.  &
       turb % model .eq. K_EPS_ZETA_F          .or.  &
       turb % model .eq. HYBRID_LES_RANS       .or.  &
       turb % model .eq. RSM_MANCEAU_HANJALIC  .or.  &
       turb % model .eq. RSM_HANJALIC_JAKIRLIC  ) then
      call Save_Scalar_Real(grid, "TurbulentKineticEnergy", plot_inside,  &
                                  turb % kin % n(-grid % n_bnd_cells),    &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentDissipation", plot_inside,    &
                                  turb % eps % n(-grid % n_bnd_cells),    &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentKineticEnergyProduction",     &
                                  plot_inside,                            &
                                  turb % p_kin(-grid % n_bnd_cells),      &
                                  f8, f9, data_offset, run)
    end if

    ! Save zeta and f22
    if(turb % model .eq. K_EPS_ZETA_F .or.  &
       turb % model .eq. HYBRID_LES_RANS) then
      v2_calc(:) = 0.0
      do c = 1, grid % n_cells
        v2_calc(c) = turb % kin % n(c) * turb % zeta % n(c)
      end do
      call Save_Scalar_Real(grid, "TurbulentQuantityV2", plot_inside,     &
                                  v2_calc (-grid % n_bnd_cells),          &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentQuantityZeta", plot_inside,   &
                                  turb % zeta % n(-grid % n_bnd_cells),   &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentQuantityF22", plot_inside,    &
                                  turb % f22  % n(-grid % n_bnd_cells),   &
                                  f8, f9, data_offset, run)
      if (heat_transfer) then
        call Save_Scalar_Real(grid, "TurbulentQuantityT2", plot_inside,   &
                                    turb % t2 % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentT2Production", plot_inside, &
                                    turb % p_t2(-grid % n_bnd_cells),     &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxX", plot_inside,    &
                                    turb % ut % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxY", plot_inside,    &
                                    turb % vt % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxZ", plot_inside,    &
                                    turb % wt % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulenQuantityAlphaL",             &
                                    plot_inside,                          &
                                    turb % alpha_l(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulenQuantityAlphaU",             &
                                    plot_inside,                          &
                                    turb % alpha_u(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
      end if
    end if

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Save_Scalar_Real(grid, "TurbulentQuantityF22", plot_inside,  &
                                  turb % f22 % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    ! Save vis and vis_t
    if(turb % model .eq. DES_SPALART .or.  &
       turb % model .eq. SPALART_ALLMARAS) then
      call Save_Scalar_Real(grid, "TurbulentViscosity", plot_inside,    &
                                  turb % vis % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "VorticityMagnitude", plot_inside,    &
                                  flow % vort(-grid % n_bnd_cells),     &
                                  f8, f9, data_offset, run)
    end if
    kin_vis_t(:) = 0.0
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      kin_vis_t   (-grid % n_bnd_cells:grid % n_cells) =  &
      turb % vis_t(-grid % n_bnd_cells:grid % n_cells) /  &
         flow % viscosity(-grid % n_bnd_cells:grid % n_cells)
      call Save_Scalar_Real(grid, "EddyOverMolecularViscosity",        &
                                  plot_inside,                         &
                                  kin_vis_t(-grid % n_bnd_cells),      &
                                  f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Save_Scalar_Real(grid, "ReynoldsStressXX", plot_inside,     &
                                  turb % uu % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressYY", plot_inside,     &
                                  turb % vv % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressZZ", plot_inside,     &
                                  turb % ww % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressXY", plot_inside,     &
                                  turb % uv % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressXZ", plot_inside,     &
                                  turb % uw % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressYZ", plot_inside,     &
                                  turb % vw % n(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
      if(heat_transfer) then
        call Save_Scalar_Real(grid, "TurbulentHeatFluxX", plot_inside,    &
                                    turb % ut % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxY", plot_inside,    &
                                    turb % vt % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxZ", plot_inside,    &
                                    turb % wt % n(-grid % n_bnd_cells),   &
                                    f8, f9, data_offset, run)
      end if
    end if

    ! Statistics for large-scale simulations of turbulence
    if(turb % statistics) then
      call Save_Vector_Real(grid, "MeanVelocity", plot_inside,          &
                                  turb % u_mean(-grid % n_bnd_cells),   &
                                  turb % v_mean(-grid % n_bnd_cells),   &
                                  turb % w_mean(-grid % n_bnd_cells),   &
                                  f8, f9, data_offset, run)
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
      call Save_Scalar_Real(grid, "MeanReynoldsStressXX", plot_inside,  &
                                  uu_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressYY", plot_inside,  &
                                  vv_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressZZ", plot_inside,  &
                                  ww_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressXY", plot_inside,  &
                                  uv_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressXZ", plot_inside,  &
                                  uw_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressYZ", plot_inside,  &
                                  vw_save(-grid % n_bnd_cells),         &
                                  f8, f9, data_offset, run)
      if(heat_transfer) then
        call Save_Scalar_Real(grid, "MeanTemperature", plot_inside,      &
                                    turb % t_mean(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
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
        call Save_Scalar_Real(grid, "MeanTurbulentQuantityT2",     &
                                    plot_inside,                   &
                                    t2_save(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxX",      &
                                    plot_inside,                   &
                                    ut_save(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxY",      &
                                    plot_inside,                   &
                                    vt_save(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxZ",      &
                                    plot_inside,                   &
                                    wt_save(-grid % n_bnd_cells),  &
                                    f8, f9, data_offset, run)
      end if

      ! Scalars
      do sc = 1, flow % n_scalars
        phi => flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        do c = 1, grid % n_cells
          phi_save(c) = turb % scalar_mean(sc, c)
        end do
        call Save_Scalar_Real(grid, name_mean, plot_inside,  &
                         phi_save(-grid % n_bnd_cells),      &
                         f8, f9, data_offset, run)
      end do
    end if

    ! Save y+ for all turbulence models
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      call Save_Scalar_Real(grid, "TurbulentQuantityYplus",            &
                                  plot_inside,                         &
                                  turb % y_plus(-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    call Save_Scalar_Real(grid, "WallDistance", plot_inside,            &
                                grid % wall_dist(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "CellDeltaMax", plot_inside,        &
                                turb % h_max(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "CellDeltaMin", plot_inside,        &
                                turb % h_min(-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "CellDeltaWall", plot_inside,       &
                                turb % h_w  (-grid % n_bnd_cells),  &
                                f8, f9, data_offset, run)

    !----------------------!
    !   Save user arrays   !
    !----------------------!
    do ua = 1, grid % n_user_arrays

      a_name = 'A_00'
      write(a_name(3:4), '(I2.2)') ua
      call Save_Scalar_Real(grid, a_name, plot_inside,                        &
                                  grid % user_array(ua,-grid % n_bnd_cells),  &
                                  f8, f9, data_offset, run)
    end do

    !----------------------!
    !                      !
    !   End of cell data   !
    !                      !
    !----------------------!
    if(run .eq. 1) then
      if(n_proc > 1 .and. this_proc .eq. 1) then
        write(f8) IN_3 // '</PCellData>' // LF
      end if
      write(f9) IN_3 // '</CellData>' // LF

      write(f9) IN_2 // '</Piece>'            // LF
      write(f9) IN_1 // '</UnstructuredGrid>' // LF

      !-------------------!
      !                   !
      !   Appended data   !
      !                   !
      !-------------------!
      write(f9) IN_0 // '<AppendedData encoding="raw">' // LF
      write(f9) '_'
    else
      write(f9) LF // IN_0 // '</AppendedData>' // LF
    end if

  end do  ! run

  write(f9) IN_0 // '</VTKFile>'          // LF
  close(f9)

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
                               extension='.vtu',  &
                               domain=domain)
      else
        call File_Mod_Set_Name(name_out_9,        &
                               processor=n,       &
                               time_step=ts,      &
                               appendix ='-bnd',  &
                               extension='.vtu',  &
                               domain=domain)
      end if
      write(f8) IN_2 // '<Piece Source="', trim(name_out_9), '"/>' // LF
    end do
    write(f8) IN_1 // '</PUnstructuredGrid>' // LF
    write(f8) IN_0 // '</VTKFile>'           // LF
    close(f8)
  end if

  call Cpu_Timer_Mod_Stop('Save_Vtu_Results')

  end subroutine
