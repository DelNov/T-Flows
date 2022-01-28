!==============================================================================!
  subroutine Save_Results(Results,  &
                          Flow, turb, Vof, swarm, ts, plot_inside, domain)
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
                      tx_save   => r_cell_18,  &
                      ty_save   => r_cell_19,  &
                      tz_save   => r_cell_20,  &
                      u_ins     => r_cell_21,  &
                      v_ins     => r_cell_22,  &
                      w_ins     => r_cell_23,  &
                      var_ins   => r_cell_24,  &
                      int_save  => i_cell_01,  &
                      type_save => i_cell_02,  &  ! cell type save array
                      offs_save => i_cell_03      ! cell offsets save array
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type)       :: Results
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: turb
  type(Vof_Type),    target :: Vof
  type(Swarm_Type),  target :: swarm
  integer                   :: ts           ! time step
  logical                   :: plot_inside  ! plot results inside?
  integer,         optional :: domain
!----------------------------------[Locals]------------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: phi
  integer(SP)              :: data_size
  integer                  :: data_offset, cell_offset, i_fac
  integer                  :: s, n, n_conns, n_polyg, sc, f8, f9, ua, run
  integer                  :: s1, s2, c1, c2, c_f, c_l
  real                     :: dist1, dist2
  character(SL)            :: name_out_8, name_out_9, name_mean, a_name
  character(SL)            :: str1, str2
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: IP           = DP  ! int. precision is double precision
  integer, parameter :: RP           = DP  ! real precision is double precision
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

  call Cpu_Timer % Start('Save_Vtu_Results')

  ! Take aliases
  Grid => Flow % pnt_grid

  if(.not. plot_inside .and. .not. Results % boundary) return

  !------------------------------------------------!
  !   Mark the beginnings and end of cell ranges   !
  !------------------------------------------------!

  ! For cells inside
  if(plot_inside) then
    c_f = 1
    c_l = Grid % n_cells
    if(.not. PLOT_BUFFERS) c_l = Grid % n_cells - Grid % Comm % n_buff_cells

  ! For boundary cells
  else
    c_f = -Grid % n_bnd_cells
    c_l = -1
    if(.not. PLOT_BUFFERS) then
      do c_f = -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_f) .eq. this_proc) exit
      end do
      do c_l = -1, -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_l) .eq. this_proc) exit
      end do
    end if
  end if

  !-------------------------------------------------------------------------!
  !   Count connections and polygons in this Grid, you will need it later   !
  !-------------------------------------------------------------------------!
  n_conns = 0
  n_polyg = 0
  if(plot_inside) then
    ! Connections
    do c1 = c_f, c_l
      n_conns = n_conns + abs(Grid % cells_n_nodes(c1))
    end do
    ! Polygons
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
        n_polyg = n_polyg + 1                   ! add one for number of polyfs
        do i_fac = 1, Grid % cells_n_faces(c1)  ! add all faces and their nodes
          s = Grid % cells_f(i_fac, c1)
          n = Grid % faces_n_nodes(s)
          n_polyg = n_polyg + 1 + n
        end do
      end if
    end do
  else
    do c2 = c_f, c_l
      n_conns = n_conns + Grid % cells_n_nodes(c2)
    end do
  end if

  call Comm_Mod_Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File % Set_Name(name_out_8,             &
                         time_step=ts,           &
                         extension='.pvtu',      &
                         domain=domain)
    call File % Set_Name(name_out_9,             &
                         processor=this_proc,    &
                         time_step=ts,           &
                         extension='.vtu',       &
                         domain=domain)
  else
    call File % Set_Name(name_out_8,             &
                         time_step=ts,           &
                         appendix ='-bnd',       &
                         extension='.pvtu',      &
                         domain=domain)
    call File % Set_Name(name_out_9,             &
                         processor=this_proc,    &
                         time_step=ts,           &
                         appendix ='-bnd',       &
                         extension='.vtu',       &
                         domain=domain)
  end if

  if(n_proc > 1 .and. this_proc .eq. 1) then
    call File % Open_For_Writing_Binary(name_out_8, f8)
  end if
  call File % Open_For_Writing_Binary(name_out_9, f9)

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

  write(str1,'(i0.0)') Grid % n_nodes
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
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3

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
      cell_offset   = cell_offset + abs(Grid % cells_n_nodes(c1))
      offs_save(c1) = cell_offset
    end do
  else
    do c2 = c_f, c_l
      cell_offset   = cell_offset + Grid % cells_n_nodes(c2)
      offs_save(c2) = cell_offset
    end do
  end if
  call Results % Save_Scalar_Int("offsets", plot_inside,  &
                                  offs_save(c_f:c_l),     &
                                  f8, f9, data_offset, 1)  ! 1 => header only

  ! Fill up an array with cell types and save the header only
  if(plot_inside) then
    do c1 = c_f, c_l
      if(Grid % cells_n_nodes(c1) .eq. 8) type_save(c1) = VTK_HEXAHEDRON
      if(Grid % cells_n_nodes(c1) .eq. 6) type_save(c1) = VTK_WEDGE
      if(Grid % cells_n_nodes(c1) .eq. 4) type_save(c1) = VTK_TETRA
      if(Grid % cells_n_nodes(c1) .eq. 5) type_save(c1) = VTK_PYRAMID
      if(Grid % cells_n_nodes(c1) .lt. 0) type_save(c1) = VTK_POLYHEDRON
    end do
  else
    do c2 = c_f, c_l
      if(Grid % cells_n_nodes(c2) .eq. 4) then
        type_save(c2) = VTK_QUAD
      else if(Grid % cells_n_nodes(c2) .eq. 3) then
        type_save(c2) = VTK_TRIANGLE
      else
        type_save(c2) = VTK_POLYGON
      end if
    end do
  end if
  call Results % Save_Scalar_Int("types", plot_inside,  &
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
    !   Save remnants of the Grid definition   !
    !------------------------------------------!
    if(run .eq. 2) then

      ! Save the nodes' coordinates
      data_size = int(Grid % n_nodes * RP * 3, SP)
      write(f9) data_size
      do n = 1, Grid % n_nodes
        write(f9) Grid % xn(n), Grid % yn(n), Grid % zn(n)
      end do

      ! Save connections
      data_size = int(n_conns * IP, SP)
      write(f9) data_size
      if(plot_inside) then
        do c1 = c_f, c_l

          ! Tetrahedral, pyramid, wedge and hexahedral cells
          if( any( Grid % cells_n_nodes(c1) .eq. (/4,5,6,8/)  ) ) then
            write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c1), c1))-1

          ! Polyhedral cells
          else if(Grid % cells_n_nodes(c1) .lt. 0) then
            write(f9) (Grid % cells_n(1:-Grid % cells_n_nodes(c1), c1))-1

          end if
        end do
      else  ! plot only boundary
        do c2 = c_f, c_l

          ! All cell types
          write(f9) (Grid % cells_n(1:Grid % cells_n_nodes(c2), c2))-1

        end do
      end if

      ! Save cell offsets
      call Results % Save_Scalar_Int("offsets", plot_inside,     &
                                      offs_save(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      ! Save cell types
      call Results % Save_Scalar_Int("types", plot_inside,       &
                                      type_save(c_f:c_l),        &
                                      f8, f9, data_offset, run)

      if(n_polyg > 0) then

        ! Write polyhedral cells' faces
        data_size = int(n_polyg * IP, SP)
        write(f9) data_size

        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            write(f9) Grid % cells_n_faces(c1)      ! write number of polygons
            do i_fac = 1, Grid % cells_n_faces(c1)  ! write nodes of each polygn
              s = Grid % cells_f(i_fac, c1)
              if(Grid % faces_s(s) .ne. 0) then  ! face has a shadow, if it ...
                s1 = s                           ! ... is closer, plot that!
                s2 = Grid % faces_s(s)
                dist1 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s1), Grid % yf(s1), Grid % zf(s1))
                dist2 = Math % Distance(                              &
                        Grid % xc(c1), Grid % yc(c1), Grid % zc(c1),  &
                        Grid % xf(s2), Grid % yf(s2), Grid % zf(s2))
                if(dist1 < dist2) s = s1
                if(dist2 < dist1) s = s2
              end if
              n = Grid % faces_n_nodes(s)
              write(f9) n, (Grid % faces_n(1:n, s))-1
            end do
          end if
        end do

        ! Write polyhedral cells' faces offsets
        data_size = int((c_l-c_f+1) * IP, SP)
        write(f9) data_size

        cell_offset = 0
        do c1 = c_f, c_l
          if(Grid % cells_n_nodes(c1) .lt. 0) then  ! found a polyhedron
            cell_offset = cell_offset + 1           ! to store number of polygs
            do i_fac = 1, Grid % cells_n_faces(c1)  ! to store all the nodes
              s = Grid % cells_f(i_fac, c1)         ! of each polygon
              n = Grid % faces_n_nodes(s)
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
      int_save(c1) = Grid % Comm % cell_proc(c1)
    end do
    do c2 = c_f, c_l
      int_save(c2) = Grid % Comm % cell_proc(c2)
    end do
    call Results % Save_Scalar_Int("Processor [1]", plot_inside,   &
                                    int_save(c_f:c_l),             &
                                    f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(c_f:c_l) = domain
      call Results % Save_Scalar_Int("Domain [1]", plot_inside,  &
                                      int_save(c_f:c_l),         &
                                      f8, f9, data_offset, run)
    end if

    !--------------!
    !   Velocity   !
    !--------------!
    call Results % Save_Vector_Real("Velocity [m/s]", plot_inside,  &
                                    Flow % u % n(c_f:c_l),          &
                                    Flow % v % n(c_f:c_l),          &
                                    Flow % w % n(c_f:c_l),          &
                                    f8, f9, data_offset, run)

    !---------------!
    !   Potential   !
    !---------------!
    call Results % Save_Scalar_Real("Potential [m^2/s]", plot_inside,  &
                                    Flow % pot % n(c_f:c_l),           &
                                    f8, f9, data_offset, run)

    !--------------------------------------!
    !   Pressure correction and pressure   !
    !--------------------------------------!
    call Results % Save_Scalar_Real("PressureCorrection [Pa]",  &
                                    plot_inside,                &
                                    Flow % pp % n(c_f:c_l),     &
                                    f8, f9, data_offset, run)
    px_save(:) = 0.0
    py_save(:) = 0.0
    pz_save(:) = 0.0
    do c1 = c_f, c_l
      px_save(c1) = Flow % pp % x(c1) * Grid % vol(c1)
      py_save(c1) = Flow % pp % y(c1) * Grid % vol(c1)
      pz_save(c1) = Flow % pp % z(c1) * Grid % vol(c1)
    end do
    call Results % Save_Vector_Real("PressureCorrectionForce [N]",  &
                                    plot_inside,                    &
                                    px_save(c_f:c_l),               &
                                    py_save(c_f:c_l),               &
                                    pz_save(c_f:c_l),               &
                                    f8, f9, data_offset, run)

    call Results % Save_Scalar_Real("Pressure [Pa]", plot_inside,  &
                                    Flow % p % n(c_f:c_l),         &
                                    f8, f9, data_offset, run)
    px_save(:) = 0.0
    py_save(:) = 0.0
    pz_save(:) = 0.0
    do c1 = c_f, c_l
      px_save(c1) = Flow % p % x(c1) * Grid % vol(c1)
      py_save(c1) = Flow % p % y(c1) * Grid % vol(c1)
      pz_save(c1) = Flow % p % z(c1) * Grid % vol(c1)
    end do
    call Results % Save_Vector_Real("PressureForce [N]", plot_inside,    &
                                    px_save(c_f:c_l),                    &
                                    py_save(c_f:c_l),                    &
                                    pz_save(c_f:c_l),                    &
                                    f8, f9, data_offset, run)

    !-----------------!
    !   Temperature   !
    !-----------------!
    if(Flow % heat_transfer) then
      call Results % Save_Scalar_Real("Temperature [K]", plot_inside,  &
                                      Flow % t % n(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      tx_save(:) = 0.0
      ty_save(:) = 0.0
      tz_save(:) = 0.0
      do c1 = c_f, c_l
        tx_save(c1) = Flow % t % x(c1)
        ty_save(c1) = Flow % t % y(c1)
        tz_save(c1) = Flow % t % z(c1)
      end do

      if(.not. Flow % mass_transfer) then
        call Flow % Grad_Variable(Flow % t)
      else
        call Vof % Calculate_Grad_Matrix_With_Front()
        call Vof % Grad_Variable_With_Front(Flow % t, Vof % t_sat)
        call Flow % Calculate_Grad_Matrix()
      end if

      call Results % Save_Vector_Real("TemperatureGradients [K/m]",  &
                                      plot_inside,                   &
                                      tx_save(c_f:c_l),              &
                                      ty_save(c_f:c_l),              &
                                      tz_save(c_f:c_l),              &
                                      f8, f9, data_offset, run)

    end if

    !-------------------------!
    !   Physical properties   !
    !-------------------------!
    call Results % Save_Scalar_Real("PhysicalDensity [kg/m^3]",      &
                                    plot_inside,                     &
                                    Flow % density(c_f:c_l),         &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("PhysicalViscosity [Pa s]",      &
                                    plot_inside,                     &
                                    Flow % viscosity(c_f:c_l),       &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("PhysicalConductivity [W/m/K]",  &
                                    plot_inside,                     &
                                    Flow % conductivity(c_f:c_l),    &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("PhysicalCapacity [J/K]",        &
                                    plot_inside,                     &
                                    Flow % capacity(c_f:c_l),        &
                                    f8, f9, data_offset, run)

    if(turb % rough_walls) then
      call Results % Save_Scalar_Real("Roughness Coefficient z_o",  &
                                      plot_inside,                  &
                                      turb % z_o_f(c_f:c_l),        &
                                      f8, f9, data_offset, run)

    end if

    !---------------------!
    !   Volume fraction   !
    !---------------------!
    if(Flow % with_interface) then
      call Results % Save_Scalar_Real("VofSharp [1]",                  &
                                      plot_inside,                     &
                                      Vof % fun % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("VofSmooth [1]",                 &
                                      plot_inside,                     &
                                      Vof % smooth % n(c_f:c_l),       &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("VofCurvature [1/m]",            &
                                      plot_inside,                     &
                                      Vof % curv(c_f:c_l),             &
                                      f8, f9, data_offset, run)
      call Results % Save_Vector_Real("VofSurfaceNormals [1]",         &
                                      plot_inside,                     &
                                      Vof % nx(c_f:c_l),               &
                                      Vof % ny(c_f:c_l),               &
                                      Vof % nz(c_f:c_l),               &
                                      f8, f9, data_offset, run)
      call Results % Save_Vector_Real("VofSurfaceTensionForce [N]",    &
                                      plot_inside,                     &
                                      Vof % surf_fx(c_f:c_l),          &
                                      Vof % surf_fy(c_f:c_l),          &
                                      Vof % surf_fz(c_f:c_l),          &
                                      f8, f9, data_offset, run)
      if (allocated(Vof % m_dot)) then
        call Results % Save_Scalar_Real("VofMassTransfer [kg/m^3/s]",  &
                                        plot_inside,                   &
                                        Vof % m_dot(c_f:c_l),          &
                                        f8, f9, data_offset, run)
      end if
    end if

    !---------------------------------------!
    !   Number of impacts and reflections   !
    !---------------------------------------!
    if(Flow % with_particles .and. .not. plot_inside) then
      call Results % Save_Scalar_Real("ParticlesReflected [1]",      &
                                      plot_inside,                   &
                                      swarm % n_reflected(c_f:c_l),  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ParticlesDeposited [1]",      &
                                      plot_inside,                   &
                                      swarm % n_deposited(c_f:c_l),  &
                                      f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      call Results % Save_Scalar_Real(phi % name, plot_inside,   &
                                      phi % n(c_f:c_l),          &
                                      f8, f9, data_offset, run)
    end do

    !-----------------!
    !   Q-criterion   !
    !-----------------!
    q_save(:) = 0.0
    do c1 = c_f, c_l
      q_save(c1) = (Flow % vort(c1)**2 - Flow % shear(c1)**2)/4.
    end do
    call Results % Save_Scalar_Real("QCriterion [1/s^2]", plot_inside,   &
                                    q_save(c_f:c_l),                     &
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
      call Results % Save_Scalar_Real("TurbulentKineticEnergy [m^2/s^2]",  &
                            plot_inside,                                   &
                            turb % kin % n(c_f:c_l),                       &
                            f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("TurbulentDissipation [m^2/s^3]",    &
                                      plot_inside,                         &
                                      turb % eps % n(c_f:c_l),             &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real(                                     &
                            "TurbulentKineticEnergyProduction [m^2/s^3]",  &
                            plot_inside,                                   &
                            turb % p_kin(c_f:c_l),                         &
                            f8, f9, data_offset, run)
    end if

    ! Save zeta and f22
    if(turb % model .eq. K_EPS_ZETA_F .or.  &
       turb % model .eq. HYBRID_LES_RANS) then
      v2_calc(:) = 0.0
      do c1 = c_f, c_l
        v2_calc(c1) = turb % kin % n(c1) * turb % zeta % n(c1)
      end do
!      call Results % Save_Scalar_Real("TurbulentQuantityV2 [m^2/s^2]",    &
!                                      plot_inside,                        &
!                                      v2_calc (c_f:c_l),                  &
!                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("TurbulentQuantityZeta [1]",        &
                                      plot_inside,                        &
                                      turb % zeta % n(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("TurbulentQuantityF22 [1]",         &
                                      plot_inside,                        &
                                      turb % f22  % n(c_f:c_l),           &
                                      f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("TurbulentQuantityT2 [K^2]",      &
                                        plot_inside,                      &
                                        turb % t2 % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("TurbulentT2Production [K^2/s]",  &
!                                        plot_inside,                      &
!                                        turb % p_t2(c_f:c_l),             &
!                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("TurbulentHeatFluxX [K m/s]",     &
                                        plot_inside,                      &
                                        turb % ut % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("TurbulentHeatFluxY [K m/s]",     &
                                        plot_inside,                      &
                                        turb % vt % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("TurbulentHeatFluxZ [K m/s]",     &
                                        plot_inside,                      &
                                        turb % wt % n(c_f:c_l),           &
                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("TurbulenQuantityAlphaL",         &
!                                        plot_inside,                      &
!                                        turb % alpha_l(c_f:c_l),          &
!                                        f8, f9, data_offset, run)
!        call Results % Save_Scalar_Real("TurbulenQuantityAlphaU",         &
!                                        plot_inside,                      &
!                                        turb % alpha_u(c_f:c_l),          &
!                                        f8, f9, data_offset, run)
      end if
    end if

    if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
      call Results % Save_Scalar_Real("TurbulentQuantityF22 [1]",  &
                                      plot_inside,                 &
                                      turb % f22 % n(c_f:c_l),     &
                                      f8, f9, data_offset, run)
    end if

    ! Save vis and vis_t
    if(turb % model .eq. DES_SPALART .or.  &
       turb % model .eq. SPALART_ALLMARAS) then
      call Results % Save_Scalar_Real("TurbulentViscosity [Pa s]",  &
                                      plot_inside,                  &
                                      turb % vis % n(c_f:c_l),      &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("VorticityMagnitude [1/s]",   &
                                      plot_inside,                  &
                                      Flow % vort(c_f:c_l),         &
                                      f8, f9, data_offset, run)
    end if

    kin_vis_t(:) = 0.0
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. HYBRID_LES_RANS     .and.  &
       turb % model .ne. DNS) then
      kin_vis_t(c_f:c_l) = turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real("EddyOverMolecularViscosity [1]",  &
                                      plot_inside,                       &
                                      kin_vis_t(c_f:c_l),                &
                                      f8, f9, data_offset, run)
    end if

    if(turb % model .eq. HYBRID_LES_RANS) then
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real("RansEddyOverMolecularViscosity [1]",  &
                                      plot_inside,                           &
                                      kin_vis_t(c_f:c_l),                    &
                                      f8, f9, data_offset, run)
      kin_vis_t(:) = 0.0
      kin_vis_t(c_f:c_l) = turb % vis_t_sgs(c_f:c_l) / Flow % viscosity(c_f:c_l)
      call Results % Save_Scalar_Real("SgsEddyOverMolecularViscosity [1]",  &
                                      plot_inside,                          &
                                      kin_vis_t(c_f:c_l),                   &
                                      f8, f9, data_offset, run)
    end if

    ! Reynolds stress models
    if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Results % Save_Scalar_Real("ReynoldsStressXX [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % uu % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ReynoldsStressYY [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % vv % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ReynoldsStressZZ [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % ww % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ReynoldsStressXY [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % uv % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ReynoldsStressXZ [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % uw % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("ReynoldsStressYZ [m^2/s^2]",  &
                                      plot_inside,                   &
                                      turb % vw % n(c_f:c_l),        &
                                      f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("TurbulentHeatFluxX [K m/s]",  &
                                        plot_inside,                   &
                                        turb % ut % n(c_f:c_l),        &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("TurbulentHeatFluxY [K m/s]",  &
                                        plot_inside,                   &
                                        turb % vt % n(c_f:c_l),        &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("TurbulentHeatFluxZ [K m/s]",  &
                                        plot_inside,                   &
                                        turb % wt % n(c_f:c_l),        &
                                        f8, f9, data_offset, run)
      end if
    end if

    ! Statistics for large-scale simulations of turbulence
    if(turb % statistics) then
      call Results % Save_Vector_Real("MeanVelocity [m/s]",      &
                                      plot_inside,               &
                                      turb % u_mean(c_f:c_l),    &
                                      turb % v_mean(c_f:c_l),    &
                                      turb % w_mean(c_f:c_l),    &
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
      call Results % Save_Scalar_Real("MeanReynoldsStressXX [m^s/s^2]",  &
                                      plot_inside,                       &
                                      uu_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("MeanReynoldsStressYY [m^s/s^2]",  &
                                      plot_inside,                       &
                                      vv_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("MeanReynoldsStressZZ [m^s/s^2]",  &
                                      plot_inside,                       &
                                      ww_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("MeanReynoldsStressXY [m^s/s^2]",  &
                                      plot_inside,                       &
                                      uv_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("MeanReynoldsStressXZ [m^s/s^2]",  &
                                      plot_inside,                       &
                                      uw_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      call Results % Save_Scalar_Real("MeanReynoldsStressYZ [m^s/s^2]",  &
                                      plot_inside,                       &
                                      vw_save(c_f:c_l),                  &
                                      f8, f9, data_offset, run)
      if(Flow % heat_transfer) then
        call Results % Save_Scalar_Real("MeanTemperature [K]",           &
                                        plot_inside,                     &
                                        turb % t_mean(c_f:c_l),          &
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
        call Results % Save_Scalar_Real("MeanTurbulentQuantityT2 [K^2]",     &
                                        plot_inside,                         &
                                        t2_save(c_f:c_l),                    &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("MeanTurbulentHeatFluxX [K m/s]",    &
                                        plot_inside,                         &
                                        ut_save(c_f:c_l),                    &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("MeanTurbulentHeatFluxY [K m/s]",    &
                                        plot_inside,                         &
                                        vt_save(c_f:c_l),                    &
                                        f8, f9, data_offset, run)
        call Results % Save_Scalar_Real("MeanTurbulentHeatFluxZ [K m/s]",    &
                                        plot_inside,                         &
                                        wt_save(c_f:c_l),                    &
                                        f8, f9, data_offset, run)
      end if

      ! Scalars
      do sc = 1, Flow % n_scalars
        phi => Flow % scalar(sc)
        name_mean = 'Mean'
        name_mean(5:8) = phi % name
        do c1 = c_f, c_l
          phi_save(c1) = turb % scalar_mean(sc, c1)
        end do
        call Results % Save_Scalar_Real(name_mean, plot_inside,    &
                                        phi_save(c_f:c_l),         &
                                        f8, f9, data_offset, run)
      end do
    end if

    ! Save y+ for all turbulence models
    if(turb % model .ne. NO_TURBULENCE_MODEL .and.  &
       turb % model .ne. DNS) then
      call Results % Save_Scalar_Real("TurbulentQuantityYplus",  &
                                      plot_inside,               &
                                      turb % y_plus(c_f:c_l),    &
                                      f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    call Results % Save_Scalar_Real("GridCellVolume [m^3]",     &
                                    plot_inside,                &
                                    Grid % vol(c_f:c_l),        &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("GridWallDistance [m]",     &
                                    plot_inside,                &
                                    Grid % wall_dist(c_f:c_l),  &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("GridCellDeltaMax [m]",     &
                                    plot_inside,                &
                                    turb % h_max(c_f:c_l),      &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("GridCellDeltaMin [m]",     &
                                    plot_inside,                &
                                    turb % h_min(c_f:c_l),      &
                                    f8, f9, data_offset, run)
    call Results % Save_Scalar_Real("GridCellDeltaWall [m]",    &
                                    plot_inside,                &
                                    turb % h_w  (c_f:c_l),      &
                                    f8, f9, data_offset, run)

    !---------------------------------------------------------------------!
    !   Variables in the first computational point, plotted at boundary   !
    !---------------------------------------------------------------------!

    ! Engage only for boundary plots (not inside means on the boundary)
    if( .not. plot_inside ) then 

      ! Initialize working variables to zero
      u_ins(:) = 0.0
      v_ins(:) = 0.0
      w_ins(:) = 0.0

      ! Copy internal values to boundary
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        if(c2 < 0) then 
          u_ins(c2) = Flow % u % n(c1)
          v_ins(c2) = Flow % v % n(c1)
          w_ins(c2) = Flow % w % n(c1)
        end if
      end do
      call Results % Save_Vector_Real("Velocity Near Wall [m/s]", plot_inside, &
                                      u_ins(c_f:c_l),                          &
                                      v_ins(c_f:c_l),                          &
                                      w_ins(c_f:c_l),                          &
                                      f8, f9, data_offset, run)

      if(turb % model .eq. K_EPS                 .or.  &
         turb % model .eq. K_EPS_ZETA_F          .or.  &
         turb % model .eq. HYBRID_LES_RANS) then

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do s = 1, Grid % n_faces
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          if(c2 < 0) then
            var_ins(c2) = Turb % kin % n(c1)
          end if
        end do

        call Results % Save_Scalar_Real("T.K.E. Near Wall [m^2/s^2]",  &
                                        plot_inside,                   &
                                        var_ins(c_f:c_l),              &
                                        f8, f9, data_offset, run)

        ! Copy internal values to boundary
        var_ins(:) = 0.0
        do s = 1, Grid % n_faces
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          if(c2 < 0) then
            var_ins(c2) = turb % y_plus(c1)
          end if
        end do

        call Results % Save_Scalar_Real("y+ Near Wall [m]",         &
                                        plot_inside,                &
                                        var_ins(c_f:c_l),           &
                                        f8, f9, data_offset, run)

        if(turb % rough_walls) then
          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = turb % id_zone(c1)
            end if
          end do

          call Results % Save_Scalar_Real("Zone [1]",                &
                                          plot_inside,               &
                                          var_ins(c_f:c_l),          &
                                          f8, f9, data_offset, run)
        end if

        do sc = 1, Flow % n_scalars
          phi => Flow % scalar(sc)
          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi % n(c1)
            end if
          end do

          call Results % Save_Scalar_Real("Scalar Near Wall",        &
                                          plot_inside,               &
                                          var_ins(c_f:c_l),          &
                                          f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi_save(c1)  ! turb % scalar_mean(sc, c1)
            end if
          end do

          call Results % Save_Scalar_Real("Mean Scalar Near Wall",  &
                                plot_inside,                        &
                                var_ins(c_f:c_l),                   &
                                f8, f9, data_offset, run)

          ! Copy internal values to boundary
          var_ins(:) = 0.0
          do s = 1, Grid % n_faces
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)
            if(c2 < 0) then
              var_ins(c2) = phi % q(c2)  ! turb % scalar_mean(sc, c1)
            end if
          end do

          call Results % Save_Scalar_Real("Wall Scalar Flux",        &
                                          plot_inside,               &
                                          var_ins(c_f:c_l),          &
                                          f8, f9, data_offset, run)

        end do
      end if
    end if

    !----------------------!
    !   Save user arrays   !
    !----------------------!
    do ua = 1, Grid % n_user_arrays

      a_name = 'A_00'
      write(a_name(3:4), '(I2.2)') ua
      call Results % Save_Scalar_Real(a_name, plot_inside,            &
                                      Grid % user_array(ua,c_f:c_l),  &
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
        call File % Set_Name(name_out_9,        &
                             processor=n,       &
                             time_step=ts,      &
                             extension='.vtu',  &
                             domain=domain)
      else
        call File % Set_Name(name_out_9,        &
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

  call Cpu_Timer % Stop('Save_Vtu_Results')

  end subroutine
