!==============================================================================!
  subroutine Results_Mod_Save(flow, turb, mult, swarm, ts, plot_inside, domain)
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
                      px_save   => r_cell_15,  &
                      py_save   => r_cell_16,  &
                      pz_save   => r_cell_17,  &
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
  integer(SP)              :: data_size
  integer                  :: data_offset, cell_offset, i_fac
  integer                  :: s, n, n_conns, n_polyg, sc, f8, f9, ua, run
  integer                  :: c1, c2, c_f, c_l
  character(SL)            :: name_out_8, name_out_9, name_mean, a_name
  character(SL)            :: str1, str2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP           = DP  ! int. precision is double precision
  integer, parameter :: RP           = DP  ! real precision is double precision
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

  call Cpu_Timer_Mod_Start('Save_Vtu_Results')

  ! Take aliases
  grid => flow % pnt_grid

  !------------------------------------------------!
  !   Mark the beginnings and end of cell ranges   !
  !------------------------------------------------!

  ! For cells inside
  if(plot_inside) then
    c_f = 1
    c_l = grid % n_cells
    if(.not. PLOT_BUFFERS) c_l = grid % n_cells - grid % comm % n_buff_cells

  ! For boundary cells
  else
    c_f = -grid % n_bnd_cells
    c_l = -1
    if(.not. PLOT_BUFFERS) then
      do c_f = -grid % n_bnd_cells, -1
        if( grid % comm % cell_proc(c_f) .eq. this_proc) exit
      end do
      do c_l = -1, -grid % n_bnd_cells, -1
        if( grid % comm % cell_proc(c_l) .eq. this_proc) exit
      end do
    end if
  end if

  !-------------------------------------------------------------------------!
  !   Count connections and polygons in this grid, you will need it later   !
  !-------------------------------------------------------------------------!
  n_conns = 0
  n_polyg = 0
  if(plot_inside) then
    ! Connections
    do c1 = c_f, c_l
      n_conns = n_conns + abs(grid % cells_n_nodes(c1))
    end do
    ! Polygons
    do c1 = c_f, c_l
      if(grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                   ! add one for number of polyfs
        do i_fac = 1, grid % cells_n_faces(c1)  ! add all faces and their nodes
          s = grid % cells_f(i_fac, c1)
          n = grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end do
  else
    do c2 = c_f, c_l
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
    write(str2,'(i0.0)') (c_l-c_f+1)
  else
    if((c_l-c_f+1) .eq. 0) then
      write(str2,'(i1)')   (c_l-c_f+1)  ! 0.0 doesn't work for zero :-/
    else
      write(str2,'(i0.0)') (c_l-c_f+1)
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
  if(plot_inside) then
    do c1 = c_f, c_l
      cell_offset   = cell_offset + abs(grid % cells_n_nodes(c1))
      offs_save(c1) = cell_offset
    end do
  else
    do c2 = c_f, c_l
      cell_offset   = cell_offset + grid % cells_n_nodes(c2)
      offs_save(c2) = cell_offset
    end do
  end if
  call Save_Scalar_Int(grid, "offsets", plot_inside,           &
                              offs_save(c_f:c_l),   &
                              f8, f9, data_offset, 1)  ! 1 => header only

  ! Fill up an array with cell types and save the header only
  if(plot_inside) then
    do c1 = c_f, c_l
      if(grid % cells_n_nodes(c1) .eq. 8) type_save(c1) = VTK_HEXAHEDRON
      if(grid % cells_n_nodes(c1) .eq. 6) type_save(c1) = VTK_WEDGE
      if(grid % cells_n_nodes(c1) .eq. 4) type_save(c1) = VTK_TETRA
      if(grid % cells_n_nodes(c1) .eq. 5) type_save(c1) = VTK_PYRAMID
      if(grid % cells_n_nodes(c1) .lt. 0) type_save(c1) = VTK_POLYHEDRON
    end do
  else
    do c2 = c_f, c_l
      if(grid % cells_n_nodes(c2) .eq. 4) then
        type_save(c2) = VTK_QUAD
      else if(grid % cells_n_nodes(c2) .eq. 3) then
        type_save(c2) = VTK_TRIANGLE
      else
        type_save(c2) = VTK_POLYGON
      end if
    end do
  end if
  call Save_Scalar_Int(grid, "types", plot_inside,             &
                              type_save(c_f:c_l),   &
                              f8, f9, data_offset, 1)  ! 1 => header only

  ! Write parts of header for polyhedral cells
  if(n_polyg > 0) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type="Int64"'        //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next


    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type="Int64"'        //  &
                      ' Name="faceoffsets"'            //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + (c_l-c_f+1) * IP  ! prepare for next

  end if

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
        do c1 = c_f, c_l

          ! Tetrahedral, pyramid, wedge and hexahedral cells
          if( any( grid % cells_n_nodes(c1) .eq. (/4,5,6,8/)  ) ) then
            write(f9) (grid % cells_n(1:grid % cells_n_nodes(c1), c1))-1

          ! Polyhedral cells
          else if(grid % cells_n_nodes(c1) .lt. 0) then
            write(f9) (grid % cells_n(1:-grid % cells_n_nodes(c1), c1))-1

          end if
        end do
      else  ! plot only boundary
        do c2 = c_f, c_l

          ! All cell types
          write(f9) (grid % cells_n(1:grid % cells_n_nodes(c2), c2))-1

        end do
      end if

      ! Save cell offsets
      call Save_Scalar_Int(grid, "offsets", plot_inside,     &
                                  offs_save(c_f:c_l),        &
                                  f8, f9, data_offset, run)
      ! Save cell types
      call Save_Scalar_Int(grid, "types", plot_inside,       &
                                  type_save(c_f:c_l),        &
                                  f8, f9, data_offset, run)

      if(n_polyg > 0) then

        ! Write polyhedral cells' faces
        data_size = n_polyg * IP
        write(f9) data_size

        do c1 = c_f, c_l
          if(grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            write(f9) grid % cells_n_faces(c1)      ! write number of polygons
            do i_fac = 1, grid % cells_n_faces(c1)  ! write nodes of each polygn
              s = grid % cells_f(i_fac, c1)
              n = grid % faces_n_nodes(s)
              write(f9) n, (grid % faces_n(1:n, s))-1
            end do
          end if
        end do

        ! Write polyhedral cells' faces offsets
        data_size = (c_l-c_f+1) * IP
        write(f9) data_size

        cell_offset = 0
        do c1 = c_f, c_l
          if(grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            cell_offset = cell_offset + 1           ! to store number of polygs
            do i_fac = 1, grid % cells_n_faces(c1)  ! to store all the nodes
              s = grid % cells_f(i_fac, c1)         ! of each polygon
              n = grid % faces_n_nodes(s)
              cell_offset = cell_offset + 1 + n
            end do
            write(f9) cell_offset

          ! Not a polyhedron, offsets are not needed and set to -1
          else
            write(f9) -1
          end if

        end do

      end if  ! n_polyg > 0

    end if

    !--------------------!
    !   Processor i.d.   !
    !--------------------!
    do c1 = c_f, c_l
      int_save(c1) = grid % comm % cell_proc(c1)
    end do
    do c2 = c_f, c_l
      int_save(c2) = grid % comm % cell_proc(c2)
    end do
    call Save_Scalar_Int(grid, "Processor", plot_inside,   &
                                int_save(c_f:c_l),         &
                                f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(c_f:c_l) = domain
      call Save_Scalar_Int(grid, "Domain", plot_inside,      &
                                  int_save(c_f:c_l),         &
                                  f8, f9, data_offset, run)
    end if

    !--------------!
    !   Velocity   !
    !--------------!
    call Save_Vector_Real(grid, "Velocity", plot_inside,   &
                                flow % u % n(c_f:c_l),     &
                                flow % v % n(c_f:c_l),     &
                                flow % w % n(c_f:c_l),     &
                                f8, f9, data_offset, run)

    !---------------!
    !   Potential   !
    !---------------!
    call Save_Scalar_Real(grid, "Potential", plot_inside,  &
                                flow % pot % n(c_f:c_l),   &
                                f8, f9, data_offset, run)

    !--------------!
    !   Pressure   !
    !--------------!
    call Save_Scalar_Real(grid, "PressureCorrection", plot_inside,   &
                                flow % pp % n(c_f:c_l),              &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "Pressure", plot_inside,             &
                                flow % p % n(c_f:c_l),               &
                                f8, f9, data_offset, run)
    px_save(:) = 0.0
    py_save(:) = 0.0
    pz_save(:) = 0.0
    do c1 = c_f, c_l
      px_save(c1) = flow % p % x(c1) * grid % vol(c1)
      py_save(c1) = flow % p % y(c1) * grid % vol(c1)
      pz_save(c1) = flow % p % z(c1) * grid % vol(c1)
    end do
    call Save_Vector_Real(grid, "PressureForce", plot_inside,  &
                                px_save(c_f:c_l),              &
                                py_save(c_f:c_l),              &
                                pz_save(c_f:c_l),              &
                                f8, f9, data_offset, run)

    !-----------------!
    !   Temperature   !
    !-----------------!
    if(heat_transfer) then
      call Save_Scalar_Real(grid, "Temperature", plot_inside,  &
                                  flow % t % n(c_f:c_l),       &
                                  f8, f9, data_offset, run)
    end if

    !-------------------------!
    !   Physical properties   !
    !-------------------------!
    call Save_Scalar_Real(grid, "PhysicalDensity", plot_inside,       &
                                flow % density(c_f:c_l),              &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "PhysicalViscosity", plot_inside,     &
                                flow % viscosity(c_f:c_l),            &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "PhysicalConductivity", plot_inside,  &
                                flow % conductivity(c_f:c_l),         &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "PhysicalCapacity", plot_inside,      &
                                flow % capacity(c_f:c_l),             &
                                f8, f9, data_offset, run)

    !---------------------!
    !   Volume fraction   !
    !---------------------!
    if(mult % model .eq. VOLUME_OF_FLUID) then
      call Save_Scalar_Real(grid, "VofSharp", plot_inside,                &
                                  mult % vof % n(c_f:c_l),                &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "VofSmooth", plot_inside,               &
                                  mult % smooth % n(c_f:c_l),             &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "VofCurvature", plot_inside,            &
                                  mult % curv(c_f:c_l),                   &
                                  f8, f9, data_offset, run)
      call Save_Vector_Real(grid, "VofSurfaceNormals", plot_inside,       &
                                  mult % nx(c_f:c_l),                     &
                                  mult % ny(c_f:c_l),                     &
                                  mult % nz(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Vector_Real(grid, "VofSurfaceTensionForce", plot_inside,  &
                                  mult % surf_fx(c_f:c_l),                &
                                  mult % surf_fy(c_f:c_l),                &
                                  mult % surf_fz(c_f:c_l),                &
                                  f8, f9, data_offset, run)
      if (allocated(mult % flux_rate)) then
        call Save_Scalar_Real(grid, "FluxRate ", plot_inside,             &
                                    mult % flux_rate(c_f:c_l),            &
                                    f8, f9, data_offset, run)
      end if
    end if

    !---------------------------------------!
    !   Number of impacts and reflections   !
    !---------------------------------------!
    if(mult % model .eq. LAGRANGIAN_PARTICLES .and. .not. plot_inside) then
      call Save_Scalar_Real(grid, "ParticlesReflected", plot_inside,  &
                                  swarm % n_reflected(c_f:c_l),       &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ParticlesDeposited", plot_inside,  &
                                  swarm % n_deposited(c_f:c_l),       &
                                  f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, flow % n_scalars
      phi => flow % scalar(sc)
      call Save_Scalar_Real(grid, phi % name, plot_inside,   &
                                  phi % n(c_f:c_l),          &
                                  f8, f9, data_offset, run)
    end do

    !-----------------!
    !   Q-criterion   !
    !-----------------!
    q_save(:) = 0.0
    do c1 = c_f, c_l
      q_save(c1) = (flow % vort(c1)**2 - flow % shear(c1)**2)/4.
    end do
    call Save_Scalar_Real(grid, "QCriterion", plot_inside,   &
                                q_save(c_f:c_l),             &
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
                                  turb % kin % n(c_f:c_l),                &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentDissipation", plot_inside,    &
                                  turb % eps % n(c_f:c_l),                &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentKineticEnergyProduction",     &
                                  plot_inside,                            &
                                  turb % p_kin(c_f:c_l),                  &
                                  f8, f9, data_offset, run)
    end if

    ! Save zeta and f22
    if(turb % model .eq. K_EPS_ZETA_F .or.  &
       turb % model .eq. HYBRID_LES_RANS) then
      v2_calc(:) = 0.0
      do c1 = c_f, c_l
        v2_calc(c1) = turb % kin % n(c1) * turb % zeta % n(c1)
      end do
      call Save_Scalar_Real(grid, "TurbulentQuantityV2", plot_inside,     &
                                  v2_calc (c_f:c_l),                      &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentQuantityZeta", plot_inside,   &
                                  turb % zeta % n(c_f:c_l),               &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "TurbulentQuantityF22", plot_inside,    &
                                  turb % f22  % n(c_f:c_l),               &
                                  f8, f9, data_offset, run)
      if (heat_transfer) then
        call Save_Scalar_Real(grid, "TurbulentQuantityT2", plot_inside,   &
                                    turb % t2 % n(c_f:c_l),               &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentT2Production", plot_inside, &
                                    turb % p_t2(c_f:c_l),                 &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxX", plot_inside,    &
                                    turb % ut % n(c_f:c_l),               &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxY", plot_inside,    &
                                    turb % vt % n(c_f:c_l),               &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxZ", plot_inside,    &
                                    turb % wt % n(c_f:c_l),               &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulenQuantityAlphaL",             &
                                    plot_inside,                          &
                                    turb % alpha_l(c_f:c_l),              &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulenQuantityAlphaU",             &
                                    plot_inside,                          &
                                    turb % alpha_u(c_f:c_l),              &
                                    f8, f9, data_offset, run)
      end if
    end if

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Save_Scalar_Real(grid, "TurbulentQuantityF22", plot_inside,  &
                                  turb % f22 % n(c_f:c_l),              &
                                  f8, f9, data_offset, run)
    end if

    ! Save vis and vis_t
    if(turb % model .eq. DES_SPALART .or.  &
       turb % model .eq. SPALART_ALLMARAS) then
      call Save_Scalar_Real(grid, "TurbulentViscosity", plot_inside,  &
                                  turb % vis % n(c_f:c_l),            &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "VorticityMagnitude", plot_inside,  &
                                  flow % vort(c_f:c_l),               &
                                  f8, f9, data_offset, run)
    end if
    kin_vis_t(:) = 0.0
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      kin_vis_t(c_f:c_l) = turb % vis_t(c_f:c_l) / flow % viscosity(c_f:c_l)
      call Save_Scalar_Real(grid, "EddyOverMolecularViscosity",  &
                                  plot_inside,                   &
                                  kin_vis_t(c_f:c_l),            &
                                  f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Save_Scalar_Real(grid, "ReynoldsStressXX", plot_inside,  &
                                  turb % uu % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressYY", plot_inside,  &
                                  turb % vv % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressZZ", plot_inside,  &
                                  turb % ww % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressXY", plot_inside,  &
                                  turb % uv % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressXZ", plot_inside,  &
                                  turb % uw % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "ReynoldsStressYZ", plot_inside,  &
                                  turb % vw % n(c_f:c_l),           &
                                  f8, f9, data_offset, run)
      if(heat_transfer) then
        call Save_Scalar_Real(grid, "TurbulentHeatFluxX", plot_inside,  &
                                    turb % ut % n(c_f:c_l),             &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxY", plot_inside,  &
                                    turb % vt % n(c_f:c_l),             &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "TurbulentHeatFluxZ", plot_inside,  &
                                    turb % wt % n(c_f:c_l),             &
                                    f8, f9, data_offset, run)
      end if
    end if

    ! Statistics for large-scale simulations of turbulence
    if(turb % statistics) then
      call Save_Vector_Real(grid, "MeanVelocity", plot_inside,  &
                                  turb % u_mean(c_f:c_l),       &
                                  turb % v_mean(c_f:c_l),       &
                                  turb % w_mean(c_f:c_l),       &
                                  f8, f9, data_offset, run)
      uu_save(:) = 0.0
      vv_save(:) = 0.0
      ww_save(:) = 0.0
      uv_save(:) = 0.0
      uw_save(:) = 0.0
      vw_save(:) = 0.0
      do c1 = c_f, c_l
        uu_save(c1) = turb % uu_res(c1) - turb % u_mean(c1) * turb % u_mean(c1)
        vv_save(c1) = turb % vv_res(c1) - turb % v_mean(c1) * turb % v_mean(c1)
        ww_save(c1) = turb % ww_res(c1) - turb % w_mean(c1) * turb % w_mean(c1)
        uv_save(c1) = turb % uv_res(c1) - turb % u_mean(c1) * turb % v_mean(c1)
        uw_save(c1) = turb % uw_res(c1) - turb % u_mean(c1) * turb % w_mean(c1)
        vw_save(c1) = turb % vw_res(c1) - turb % v_mean(c1) * turb % w_mean(c1)
      end do
      call Save_Scalar_Real(grid, "MeanReynoldsStressXX", plot_inside,  &
                                  uu_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressYY", plot_inside,  &
                                  vv_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressZZ", plot_inside,  &
                                  ww_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressXY", plot_inside,  &
                                  uv_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressXZ", plot_inside,  &
                                  uw_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      call Save_Scalar_Real(grid, "MeanReynoldsStressYZ", plot_inside,  &
                                  vw_save(c_f:c_l),                     &
                                  f8, f9, data_offset, run)
      if(heat_transfer) then
        call Save_Scalar_Real(grid, "MeanTemperature", plot_inside,     &
                                    turb % t_mean(c_f:c_l),             &
                                    f8, f9, data_offset, run)
        t2_save(:) = 0.0
        ut_save(:) = 0.0
        vt_save(:) = 0.0
        wt_save(:) = 0.0
        do c1 = c_f, c_l
          t2_save(c1) = turb % t2_res(c1) - turb % t_mean(c1)*turb % t_mean(c1)
          ut_save(c1) = turb % ut_res(c1) - turb % u_mean(c1)*turb % t_mean(c1)
          vt_save(c1) = turb % vt_res(c1) - turb % v_mean(c1)*turb % t_mean(c1)
          wt_save(c1) = turb % wt_res(c1) - turb % w_mean(c1)*turb % t_mean(c1)
        end do
        call Save_Scalar_Real(grid, "MeanTurbulentQuantityT2",     &
                                    plot_inside,                   &
                                    t2_save(c_f:c_l),              &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxX",      &
                                    plot_inside,                   &
                                    ut_save(c_f:c_l),              &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxY",      &
                                    plot_inside,                   &
                                    vt_save(c_f:c_l),              &
                                    f8, f9, data_offset, run)
        call Save_Scalar_Real(grid, "MeanTurbulentHeatFluxZ",      &
                                    plot_inside,                   &
                                    wt_save(c_f:c_l),              &
                                    f8, f9, data_offset, run)
      end if

      ! Scalars
      do sc = 1, flow % n_scalars
        phi => flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        do c1 = c_f, c_l
          phi_save(c1) = turb % scalar_mean(sc, c1)
        end do
        call Save_Scalar_Real(grid, name_mean, plot_inside,  &
                         phi_save(c_f:c_l),                  &
                         f8, f9, data_offset, run)
      end do
    end if

    ! Save y+ for all turbulence models
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      call Save_Scalar_Real(grid, "TurbulentQuantityYplus",  &
                                  plot_inside,               &
                                  turb % y_plus(c_f:c_l),    &
                                  f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    call Save_Scalar_Real(grid, "GridWallDistance", plot_inside,   &
                                grid % wall_dist(c_f:c_l),         &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "GridCellDeltaMax", plot_inside,   &
                                turb % h_max(c_f:c_l),             &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "GridCellDeltaMin", plot_inside,   &
                                turb % h_min(c_f:c_l),             &
                                f8, f9, data_offset, run)
    call Save_Scalar_Real(grid, "GridCellDeltaWall", plot_inside,  &
                                turb % h_w  (c_f:c_l),             &
                                f8, f9, data_offset, run)

    !----------------------!
    !   Save user arrays   !
    !----------------------!
    do ua = 1, grid % n_user_arrays

      a_name = 'A_00'
      write(a_name(3:4), '(I2.2)') ua
      call Save_Scalar_Real(grid, a_name, plot_inside,            &
                                  grid % user_array(ua,c_f:c_l),  &
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
