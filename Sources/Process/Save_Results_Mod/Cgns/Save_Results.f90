!==============================================================================!
  subroutine Save_Results(flow, turb, name_save)
!------------------------------------------------------------------------------!
!   Creates save file and adds fields to existing grid cgns                    !
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
  character(len=80)        :: store_name, name_out
  integer                  :: base
  integer                  :: block
  integer                  :: solution
  integer                  :: field
  integer                  :: c
!==============================================================================!

  call Cpu_Timer_Mod_Start('Save_Cgns_Results')

  ! Take aliases
  grid => flow % pnt_grid

  ! Store the name
  store_name = problem_name

  problem_name = name_save

  call Name_File(0, name_out, ".cgns")

  if (this_proc .lt. 2) print *, "# subroutine Save_Cgns_Results"

  !-------------------------------------------!
  !   Write a mesh (if not already written)   !
  !-------------------------------------------!

  call Save_Cgns_Cells(grid, this_proc)

  !--------------------------!
  !   Open file for modify   !
  !--------------------------!
  file_mode = CG_MODE_MODIFY
  call Cgns_Mod_Open_File(name_out, file_mode)


  call Cgns_Mod_Initialize_Counters

  ! Count number of 3d cell type elements
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    if(grid % cells_n_nodes(c) .eq. 8) cnt_hex = cnt_hex + 1
    if(grid % cells_n_nodes(c) .eq. 6) cnt_wed = cnt_wed + 1
    if(grid % cells_n_nodes(c) .eq. 5) cnt_pyr = cnt_pyr + 1
    if(grid % cells_n_nodes(c) .eq. 4) cnt_tet = cnt_tet + 1
  end do

  !-----------------!
  !                 !
  !   Bases block   !
  !                 !
  !-----------------!
  n_bases = 1
  allocate(cgns_base(n_bases))

  base = 1
  cgns_base(base) % name = "Base 1"
  cgns_base(base) % cell_dim = 3
  cgns_base(base) % phys_dim = 3

  !-----------------!
  !                 !
  !   Zones block   !
  !                 !
  !-----------------!

  cgns_base(base) % n_blocks = 1
  allocate(cgns_base(base) % block(cgns_base(base) % n_blocks))

  block = 1
  cgns_base(base) % block(block) % name = "Zone 1"
  cgns_base(base) % block(block) % mesh_info(1) = grid % n_nodes
  cgns_base(base) % block(block) % mesh_info(2) = grid % n_cells - &
                                                  grid % comm % n_buff_cells
  cgns_base(base) % block(block) % mesh_info(3) = 0

  !--------------------!
  !                    !
  !   Solution block   !
  !                    !
  !--------------------!

  cgns_base(base) % block(block) % n_solutions = 1

  allocate(cgns_base(base) % block(block) % solution( &
    cgns_base(base) % block(block) % n_solutions))
  solution = 1

  cgns_base(base) % block(block) % solution(solution) % name = "FlowSolution"
  cgns_base(base) % block(block) % solution(solution) % sol_type = CellCenter

  call Cgns_Mod_Write_Solution_Info(base, block, solution)

  !-----------------!
  !                 !
  !   Field block   !
  !                 !
  !-----------------!

  !---------------------------------------------!
  !   Copied code below from Save_Vtu_Results   !
  !---------------------------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                            flow % u % n(1), "VelocityX")
  call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                            flow % v % n(1), "VelocityY")
  call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                            flow % w % n(1), "VelocityZ")
  !--------------!
  !   Pressure   !
  !--------------!
  call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                            flow % p % n(1), "Pressure")
  !-----------------!
  !   Temperature   !
  !-----------------!
  if(heat_transfer) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              flow % t % n(1), "Temperature")
  end if

  !--------------------------!
  !   Turbulent quantities   !
  !--------------------------!

  ! Kin and Eps and T2 for K_EPS (still a little confusing)
  if(turbulence_model .eq. K_EPS                 .or.  &
     turbulence_model .eq. K_EPS_ZETA_F          .or.  &
     turbulence_model .eq. HYBRID_LES_RANS       .or.  &
     turbulence_model .eq. RSM_MANCEAU_HANJALIC  .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC       ) then

    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % kin % n(1), "TurbulentKineticEnergy")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % eps % n(1), "TurbulentDissipation")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              p_kin(1), "TurbulentKineticEnergyProduction")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % alpha1(1), "TurbulentQauntityAlpha1")

    if(heat_transfer .and. turbulence_model .eq. K_EPS) then
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                turb % t2 % n(1),  "TurbulentQuantityT2")
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                p_t2(1),    "TurbulentT2Production")
    end if
   
  end if

  ! Zeta, v2 and f22 and T2
  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      v2_calc(c) = turb % kin % n(c) * turb % zeta % n(c)
    end do
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              v2_calc(1),  "TurbulentQuantityV2")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % zeta % n(1), "TurbulentQuantityZeta")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % f22 % n(1),  "TurbulentQuantityF22")
 
    if(heat_transfer) then
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                turb % t2 % n(1),  "TurbulentQuantityT2")
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                p_t2(1),    "TurbulentT2Production")
    end if
    
  end if

  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % f22 % n(1),  "TurbulentQuantityF22")
  end if

  ! Vis and Turbulent Vicosity_t
  if(turbulence_model .eq. DES_SPALART .or.  &
     turbulence_model .eq. SPALART_ALLMARAS) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % vis % n(1),"TurbulentViscosity")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              flow % vort(1),"VorticityMagnitude")
  end if
  if(turbulence_model .ne. NONE) then
    kin_vis_t(1) = vis_t(1) / viscosity
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              kin_vis_t(1),"EddyOverMolecularViscosity")
  end if

  ! Reynolds stress models
  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % uu % n(1),"ReynoldsStressXX")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % vv % n(1),"ReynoldsStressYY")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % ww % n(1),"ReynoldsStressZZ")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % uv % n(1),"ReynoldsStressXY")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % uw % n(1),"ReynoldsStressXZ")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % vw % n(1),"ReynoldsStressYZ")
  end if

  ! Statistics for large-scale simulations of turbulence
  if(turbulence_model .eq. LES_SMAGORINSKY    .or.  &
     turbulence_model .eq. LES_DYNAMIC        .or.  &
     turbulence_model .eq. LES_WALE           .or.  &
     turbulence_model .eq. DNS                .or.  &
     turbulence_model .eq. DES_SPALART        .or.  &
     turbulence_model .eq. HYBRID_LES_PRANDTL .or.  &
     turbulence_model .eq. HYBRID_LES_RANS) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % u_mean(1),"MeanVelocityX")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % v_mean(1),"MeanVelocityY")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % w_mean(1),"MeanVelocityZ")
    do c = 1, grid % n_cells
      uu_save(c) = turb % uu_res(c) - turb % u_mean(c) * turb % u_mean(c)
      vv_save(c) = turb % vv_res(c) - turb % v_mean(c) * turb % v_mean(c)
      ww_save(c) = turb % ww_res(c) - turb % w_mean(c) * turb % w_mean(c)
      uv_save(c) = turb % uv_res(c) - turb % u_mean(c) * turb % v_mean(c)
      uw_save(c) = turb % uw_res(c) - turb % u_mean(c) * turb % w_mean(c)
      vw_save(c) = turb % vw_res(c) - turb % v_mean(c) * turb % w_mean(c)
    end do

    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              uu_save(1),"ReynoldsStressXX")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              vv_save(1),"ReynoldsStressYY")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              ww_save(1),"ReynoldsStressZZ")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              uv_save(1),"ReynoldsStressXY")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              uw_save(1),"ReynoldsStressXZ")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              vw_save(1),"ReynoldsStressYZ")
    if(heat_transfer) then
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                               turb % t_mean(1), "TemperatureMean")
      do c = 1, grid % n_cells
        t2_save(c) = turb % t2_res(c) - turb % t_mean(c) * turb % t_mean(c)
        ut_save(c) = turb % ut_res(c) - turb % u_mean(c) * turb % t_mean(c)
        vt_save(c) = turb % vt_res(c) - turb % v_mean(c) * turb % t_mean(c)
        wt_save(c) = turb % wt_res(c) - turb % w_mean(c) * turb % t_mean(c)
      end do

      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                t2_save(1),"TemperatureFluctuations")
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                ut_save(1),"TurbulentHeatFluxX")
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                vt_save(1),"TurbulentHeatFluxY")
      call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                                wt_save(1),"TurbulentHeatFluxZ")
    end if
  end if

  ! Save y+ for all turbulence models
  if(turbulence_model .ne. NONE) then
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              y_plus(1),"TurbulentQuantityYplus")
  end if

  ! Wall distance and delta
  if (.not. permanent_fields_written) then ! [actual write]

    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              grid % wall_dist(1),"WallDistance")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % h_max(1),"CellDeltaMax")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % h_min(1),"CellDeltaMin")
    call Cgns_Mod_Write_Field(base, block, solution, field, grid, &
                              turb % h_w(1),"CellDeltaWall")
    permanent_fields_written = .true.
  else ! [link]
    call Write_Link_To_Field(base, block, solution, "WallDistance")
    call Write_Link_To_Field(base, block, solution, "CellDeltaMax")
    call Write_Link_To_Field(base, block, solution, "CellDeltaMin")
    call Write_Link_To_Field(base, block, solution, "CellDeltaWall")
  end if

  !----------------------!
  !   Save user arrays   !
  !----------------------!
  ! call User_Mod_Save_Cgns_Results(base, block, solution, field, grid)

  !----------------------------!
  !   Add info on dimensions   !
  !----------------------------!
  call Write_Dimensions_Info(base, block)

  ! Close DB
  call Cgns_Mod_Close_File

  if (this_proc < 2) &
    print *, "# Added fields to ", trim(name_out)

  deallocate(cgns_base)

  ! Restore the name
  problem_name = store_name

  call Cpu_Timer_Mod_Stop('Save_Cgns_Results')

  end subroutine
