!==============================================================================!
  subroutine Save_Vtu_Fields(Results, Turb, Flow, Grid,  &
                             plot_inside, domain)
!------------------------------------------------------------------------------!
!>  The Save_Vtu_Fields subroutine is a comprehensive and central
!>  function for writing simulation results in the VTU (VTK Unstructured
!>  Grid) file format, which is compatible with visualization tools like
!>  VisIt and ParaView. This subroutine handles the intricacies of
!>  formatting and outputting a wide range of simulation data, making it
!>  crucial for post-processing and analysis
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Flexible output: Capable of outputting results for both internal         !
!     domain cells (plot_inside) and boundary cells. This flexibility allows   !
!     users to focus on specific areas of interest in the simulation           !
!   * Variable handling: Manages a variety of simulation variables             !
!     including velocity, pressure, temperature, turbulent quantities, and     !
!     more. It ensures these variables are properly formatted and written to   !
!     the VTU file.                                                            !
!   * Grid and cell data management: Calculates and prepares grid and cell     !
!     data, including coordinates, connectivity, offsets, and types, which     !
!     are essential components of the VTU format.                              !
!   * Polyhedral cell support: Specifically handles polyhedral cells, a        !
!     complex aspect often encountered in unstructured grid simulations.       !
!   * Scalar and vector data writing: Incorporates functions to write both     !
!     scalar and vector data, catering to different types of physical          !
!     quantities in the simulation.                                            !
!   * Integration with turbulence modeling: Seamlessly integrates with         !
!     various turbulence models in T-Flows, ensuring that relevant turbulent   !
!     quantities are included in the results.                                  !
!   * Boundary condition processing: Includes special handling for boundary    !
!     cells, enabling detailed analysis of boundary conditions and near-wall   !
!     phenomena.                                                               !
!   * Support for large-scale turbulence simulation statistics: Facilitates    !
!     the output of statistical data crucial for understanding turbulence,     !
!     particularly in large-scale simulations.                                 !
!   * Python script generation: Generates Python scripts for extracting        !
!     boundary conditions, enhancing post-processing capabilities with         !
!     external tools.                                                          !
!   * File management: Handles the creation and management of VTU and PVTU     !
!     files, ensuring proper file structure and format for visualization       !
!     tools.                                                                   !
!------------------------------------------------------------------------------!
!   Workflow                                                                   !
!                                                                              !
!   * Initialization: Sets up precision, aliases, grid pointers, and checks    !
!     conditions for boundary plotting.                                        !
!   * Buffer allocation and grid data preparation: Allocates memory for        !
!     various buffers and prepares grid-related data, including node           !
!     coordinates and cell connections.                                        !
!   * File creation and header writing: Creates VTU and PVTU files, writing    !
!     necessary headers for unstructured grid data.                            !
!   * Data writing in two sweeps:                                              !
!     - First sweep: Writes header information for various data types          !
!       including nodes, cells, and results.                                   !
!     - Second sweep: Appends actual data to the files, including              !
!       coordinates, connectivity, offsets, types, and simulation results.     !
!   * Results output: Outputs a wide range of simulation results, including    !
!     physical properties, turbulence data, and specific cell data.            !
!   * Polyhedral cells handling: Manages data specific to polyhedral cells,    !
!     if present.                                                              !
!   * Footer writing and file closure: Finalizes the VTU and PVTU files by     !
!     writing footers and closing files.                                       !
!   * Python script generation for boundary conditions: Optionally             !
!     generates Python scripts for extracting and analyzing boundary           !
!     conditions.                                                              !
!   * Cleanup: Disconnects from various data buffers and stops the profiler.   !
!------------------------------------------------------------------------------!
  implicit none
!--------------------------------[Arguments]-----------------------------------!
  class(Results_Type), target :: Results      !! parent class
  type(Turb_Type),     target :: Turb         !! turbulence models
  type(Field_Type),    target :: Flow         !! flow field
  type(Grid_Type)             :: Grid         !! computational grid
  logical                     :: plot_inside  !! true to plots inside,
                                              !! false to plot on the boundary
  integer,           optional :: domain       !! computational domain
!----------------------------------[Locals]------------------------------------!
  type(Var_Type),  pointer     :: phi
  integer(SP)                  :: data_size
  integer                      :: data_offset, cell_offset, i_fac, reg
  integer                      :: s, n, n_conns, n_polyg, sc, f8, f9, run
  integer                      :: s1, s2, c1, c2, c_f, c_l
  real                         :: dist1, dist2
  character(SL)                :: name_out_8, name_out_9, name_mean
  character(SL)                :: str1, str2, str_time, str_var
  integer, pointer, contiguous :: int_save(:), type_save(:), offs_save(:)
  real,    pointer, contiguous :: save_01(:), save_02(:), save_03(:)
  real,    pointer, contiguous :: save_04(:), save_05(:), save_06(:)
  real,    pointer, contiguous :: var_ins(:)
  real,    pointer, contiguous ::             kin_vis_t(:)
  logical, pointer             :: units
!------------------------------[Local parameters]------------------------------!
  logical, parameter :: PLOT_BUFFERS = .false.  ! .true. is good for debugging
!==============================================================================!

  call Profiler % Start('Save_Vtu_Results')

  ! Set precision for plotting (intp and floatp variables)
  call Vtk_Mod_Set_Precision()

  ! Take aliases
  units => Results % units

  if(.not. plot_inside .and. .not. Results % boundary) return

  call Work % Connect_Int_Cell(int_save, type_save, offs_save)
  call Work % Connect_Real_Cell(save_01, save_02, save_03,  &
                                save_04, save_05, save_06)
  call Work % Connect_Real_Cell(                  kin_vis_t          )

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
        if( Grid % Comm % cell_proc(c_f) .eq. This_Proc()) exit
      end do
      do c_l = -1, -Grid % n_bnd_cells, -1
        if( Grid % Comm % cell_proc(c_l) .eq. This_Proc()) exit
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

  call Global % Wait

  !--------------------------------------!
  !                                      !
  !   Create .pvtu file and .vtu files   !
  !                                      !
  !--------------------------------------!
  if(plot_inside) then
    call File % Set_Name(name_out_8,                    &
                         time_step = Time % Curr_Dt(),  &
                         extension = '.pvtu',           &
                         domain    = domain)
    call File % Set_Name(name_out_9,                              &
                         processor = (/This_Proc(), N_Procs()/),  &
                         time_step = Time % Curr_Dt(),            &
                         extension = '.vtu',                      &
                         domain    = domain)
  else
    call File % Set_Name(name_out_8,                    &
                         time_step = Time % Curr_Dt(),  &
                         appendix  = '-bnd',            &
                         extension = '.pvtu',           &
                         domain=domain)
    call File % Set_Name(name_out_9,                              &
                         processor = (/This_Proc(), N_Procs()/),  &
                         time_step = Time % Curr_Dt(),            &
                         appendix  = '-bnd',                      &
                         extension = '.vtu',                      &
                         domain    = domain)
  end if

  if(Parallel_Run() .and. First_Proc()) then
    call File % Open_For_Writing_Binary(name_out_8, f8)
  end if
  call File % Open_For_Writing_Binary(name_out_9, f9)

  write(str_time,'(E16.9)')Time % Get_Time()
  !------------!
  !            !
  !   Header   !
  !            !
  !------------!
  if(Parallel_Run() .and. First_Proc())  then
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
  if(Parallel_Run() .and. First_Proc())  then
    write(f8) IN_3 // '<PPoints>' // LF
    write(f8) IN_4 // '<PDataArray type='//floatp  //  &
                      ' NumberOfComponents="3"/>'  // LF
    write(f8) IN_3 // '</PPoints>' // LF
  end if

  write(str1, '(i1)') data_offset
  write(f9) IN_3 // '<Points>'                        // LF
  write(f9) IN_4 // '<DataArray type='//floatp        //  &
                    ' NumberOfComponents="3"'         //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
  write(f9) IN_4 // '</DataArray>' // LF
  write(f9) IN_3 // '</Points>'    // LF
  data_offset = data_offset + SP + Grid % n_nodes * RP * 3

  !-----------!
  !   Cells   !
  !-----------!
  write(f9) IN_3 // '<Cells>' // LF

  ! Cells' nodes
  write(str1, '(i0.0)') data_offset
  write(f9) IN_4 // '<DataArray type='//intp          //  &
                    ' Name="connectivity"'            //  &
                    ' format="appended"'              //  &
                    ' offset="' // trim(str1) //'">'  // LF
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
  call Results % Save_Vtu_Scalar_Int("offsets", plot_inside,  &
                                     offs_save(c_f:c_l),      &
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
  call Results % Save_Vtu_Scalar_Int("types", plot_inside,  &
                                     type_save(c_f:c_l),    &
                                     f8, f9, data_offset, 1)  ! 1 => header only

  ! Write parts of header for polyhedral cells
  if(n_polyg > 0) then

    ! Write polyhedral cells' faces
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
                      ' Name="faces"'                  //  &
                      ' format="appended"'             //  &
                      ' offset="' // trim(str1) //'">' // LF
    write(f9) IN_4 // '</DataArray>' // LF
    data_offset = data_offset + SP + n_polyg * IP  ! prepare for next


    ! Write polyhedral cells' faces offsets
    write(str1, '(i0.0)') data_offset
    write(f9) IN_4 // '<DataArray type='//intp         //  &
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
  if(Parallel_Run() .and. First_Proc())  then
    write(f8) IN_3 // '<PCellData>' // LF
  end if
  write(f9) IN_3 // '<CellData>' // LF

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
      call Results % Save_Vtu_Scalar_Int("offsets", plot_inside,    &
                                         offs_save(c_f:c_l),        &
                                         f8, f9, data_offset, run)
      ! Save cell types
      call Results % Save_Vtu_Scalar_Int("types", plot_inside,      &
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
    !   Cell processor   !
    !--------------------!
    do c1 = c_f, c_l
      int_save(c1) = Grid % Comm % cell_proc(c1)
    end do
    do c2 = c_f, c_l
      int_save(c2) = Grid % Comm % cell_proc(c2)
    end do

    str_var = Results % Var_Name("Grid Processor","[1]", units)
    call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                    int_save(c_f:c_l),               &
                                    f8, f9, data_offset, run)

    !-----------------!
    !   Cell thread   !
    !-----------------!
    int_save(c_f:c_l) = Grid % Omp % cell_thread(c_f:c_l)
    str_var = Results % Var_Name("Grid Thread","[1]", units)
    call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                       int_save(c_f:c_l),            &
                                       f8, f9, data_offset, run)

    !-------------------!
    !   Domain number   !
    !-------------------!
    if(present(domain)) then
      int_save(c_f:c_l) = domain
      str_var = Results % Var_Name("Grid Domain","[1]", units)
      call Results % Save_Vtu_Scalar_Int(trim(str_var), plot_inside,   &
                                         int_save(c_f:c_l),            &
                                         f8, f9, data_offset, run)
    end if

    !----------------------!
    !   Boundary regions   !
    !----------------------!
    if(.not. plot_inside) then
      int_save(c_f:c_l) = Grid % region % at_cell(c_f:c_l)
      str_var = Results % Var_Name("Boundary Condition","[1]", units)
      call Results % Save_Vtu_Scalar_Int(trim(str_var),plot_inside,    &
                                         int_save(c_f:c_l),            &
                                         f8, f9, data_offset, run)
    end if

    !--------------!
    !   Velocity   !
    !--------------!
    str_var = Results % Var_Name("Velocity","[m/s]", units)
    call Results % Save_Vtu_Vector_Real(trim(str_var), plot_inside,   &
                                        Flow % u % n(c_f:c_l),        &
                                        Flow % v % n(c_f:c_l),        &
                                        Flow % w % n(c_f:c_l),        &
                                        f8, f9, data_offset, run)

    !--------------------------------------!
    !   Pressure correction and pressure   !
    !--------------------------------------!
    str_var = Results % Var_Name("Pressure Correction","[Pa]", units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),            &
                                        plot_inside,              &
                                        Flow % pp % n(c_f:c_l),   &
                                        f8, f9, data_offset, run)

    str_var = Results % Var_Name("Pressure","[Pa]", units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,  &
                                        Flow % p % n(c_f:c_l),         &
                                        f8, f9, data_offset, run)
    !-----------------!
    !   Temperature   !
    !-----------------!
    if(Flow % heat_transfer) then
      str_var = Results % Var_Name("Temperature","[K]", units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Flow % t % n(c_f:c_l),        &
                                          f8, f9, data_offset, run)
    end if

    !-------------------------!
    !   Physical properties   !
    !-------------------------!
    str_var = Results % Var_Name("Physical Density","[kg/m^3]", units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % density(c_f:c_l),        &
                                        f8, f9, data_offset, run)
    str_var = Results % Var_Name("Physical Viscosity","[Pa s]", units)
    call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                        plot_inside,                    &
                                        Flow % viscosity(c_f:c_l),      &
                                        f8, f9, data_offset, run)
    if(Flow % heat_transfer) then
      str_var = Results % Var_Name("Physical Conductivity","[W/m/K]", units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                          plot_inside,                    &
                                          Flow % conductivity(c_f:c_l),   &
                                          f8, f9, data_offset, run)
      str_var = Results % Var_Name("Physical Capacity","[J/kg/K]", units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var),                  &
                                          plot_inside,                    &
                                          Flow % capacity(c_f:c_l),       &
                                          f8, f9, data_offset, run)
    end if

    !------------------!
    !   Save scalars   !
    !------------------!
    do sc = 1, Flow % n_scalars
      phi => Flow % scalar(sc)
      call Results % Save_Vtu_Scalar_Real(phi % name, plot_inside,   &
                                          phi % n(c_f:c_l),          &
                                          f8, f9, data_offset, run)
    end do

    !--------------------------!
    !   Turbulent quantities   !
    !--------------------------!

   kin_vis_t(:) = 0.0
    if(Turb % model .ne. NO_TURBULENCE_MODEL) then
      kin_vis_t(c_f:c_l) = Turb % vis_t(c_f:c_l) / Flow % viscosity(c_f:c_l)
      str_var = Results % Var_Name("Eddy Over Molecular Viscosity","[1]",  &
                                   units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,      &
                                          kin_vis_t(c_f:c_l),              &
                                          f8, f9, data_offset, run)
    end if

    ! Wall distance and delta, important for all models
    if(Turb % model .ne. NO_TURBULENCE_MODEL) then
      str_var = Results % Var_Name("Grid Wall Distance","[m]", units)
      call Results % Save_Vtu_Scalar_Real(trim(str_var), plot_inside,   &
                                          Grid % wall_dist(c_f:c_l),    &
                                          f8, f9, data_offset, run)
    end if

    !----------------------!
    !                      !
    !   End of cell data   !
    !                      !
    !----------------------!
    if(run .eq. 1) then
      if(Parallel_Run() .and. First_Proc()) then
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
  if(Parallel_Run() .and. First_Proc()) then
    do n = 1, N_Procs()
      if(plot_inside) then
        call File % Set_Name(name_out_9,                    &
                             processor = (/n, N_Procs()/),  &
                             time_step = Time % Curr_Dt(),  &
                             extension = '.vtu',            &
                             domain    = domain)
      else
        call File % Set_Name(name_out_9,                    &
                             processor = (/n, N_Procs()/),  &
                             time_step = Time % Curr_Dt(),  &
                             appendix  = '-bnd',            &
                             extension = '.vtu',            &
                             domain    =  domain)
      end if
      write(f8) IN_2 // '<Piece Source="', trim(name_out_9), '"/>' // LF
    end do
    write(f8) IN_1 // '</PUnstructuredGrid>' // LF
    write(f8) IN_0 // '</VTKFile>'           // LF
    close(f8)
  end if

  call Work % Disconnect_Int_Cell(int_save, type_save, offs_save)
  call Work % Disconnect_Real_Cell(save_01, save_02, save_03,  &
                                   save_04, save_05, save_06)
  call Work % Disconnect_Real_Cell(                  kin_vis_t          )

  call Profiler % Stop('Save_Vtu_Results')

  end subroutine
