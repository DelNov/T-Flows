!==============================================================================!
  subroutine Form_Momentum_Matrix(Process, Grid, Flow, Turb,  &
                                  visc_eff, urf, dt)
!------------------------------------------------------------------------------!
!   Momentum matrix is formed in the following steps:
!
!   * Physical properties setup
!     - An alias for density is defined
!     - Effective viscosity is computed as the sum of laminar and turbulent
!   * Matrix is initialized to zero
!   * Matrix coefficients are computed
!     - Viscous coefficients inside the domain first
!     - Upwind blending coefficients in the domain follow
!     - Viscous coefficients on the boundary
!     - Upwind blending coefficients on the boundary
!   * Diagonal matrix entry for the unsteady term is formed next
!   * Entries for pressure matrix are stored
!   * Matrix is under-relaxed
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)                 :: Process
  type(Grid_Type), intent(in), target :: Grid
  type(Field_Type),            target :: Flow
  type(Turb_Type),             target :: Turb
  real                                :: visc_eff(-Grid % n_bnd_cells &
                                                  :Grid % n_cells)
  real                                :: urf
  real,  optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  integer,   contiguous, pointer :: row(:)
  real,      contiguous, pointer :: dens(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12, a21, fl, cfs, w1, w2
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Momentum_Matrix')

  !-----------------------!
  !   Take some aliases   !
  !-----------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % A % dia
  pos => Flow % Nat % A % pos
  row => Flow % Nat % A % row
  fc  => Flow % Nat % A % fc
  nz  =  Flow % Nat % A % nonzeros

  Assert(urf > 0.0)

  !-------------------------------!
  !                               !
  !   Physical properties setup   !
  !                               !
  !-------------------------------!

  !------------------------!
  !   Initialize density   !
  !------------------------!
  dens => Flow % density

  !-----------------------------------------------------------!
  !   Start by copying molecular viscosity to the effective   !
  !-----------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    visc_eff(c) = Flow % viscosity(c)
  end do
  !$tf-acc loop end

  !-------------------------------------------------------------!
  !   If there is a turbulence model, add turbulent viscosity   !
  !-------------------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then

    ! Inside the domain, add turbulent viscosity
    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      visc_eff(c) = visc_eff(c) + Turb % vis_t(c)
    end do
    !$tf-acc loop end

    ! On the walls, add wall viscosity and on other
    ! boundaries just add turbulent viscosity
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL    .or.  &
         Grid % region % type(reg) .eq. WALLFL) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c1 = Grid % faces_c(1,s)   ! inside cell
          visc_eff(c1) = Turb % vis_w(c1)
        end do
        !$tf-acc loop end

      end if
    end do

  end if

  !---------------------------------------!
  !                                       !
  !   Initialize matrix entries to zero   !
  !                                       !
  !---------------------------------------!

  !$tf-acc loop begin
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$tf-acc loop end

  !---------------------------------------!
  !                                       !
  !   Compute neighbouring coefficients   !
  !                                       !
  !---------------------------------------!

  !----------------------------------------------!
  !   Viscosity coefficients inside the domain   !
  !----------------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = Grid % cells_i_cells(c1),  &  ! first inside cell
               Grid % cells_n_cells(c1)

      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      w1 = Grid % f(s)
      if(c1.gt.c2) w1 = 1.0 - w1
      w2 = 1.0 - w1

      a12 = (w1 * visc_eff(c1) + w2 * visc_eff(c2)) * fc(s)
      a21 = a12

      if(c1 .lt. c2) then
        val(pos(1,s)) = -a12
        val(pos(2,s)) = -a21
      end if

      ! Update only diaginal at c1 to avoid race conditions
      val(dia(c1)) = val(dia(c1)) + a12

    end do

  end do
  !$tf-acc loop end

  !---------------------------------------!
  !   Upwind blending inside the domain   !
  !---------------------------------------!
  if(Flow % u % blend_matrix) then

    !$tf-acc loop begin
    do c1 = Cells_In_Domain()  ! all present

      do i_cel = Grid % cells_i_cells(c1),  &  ! first inside neighbour
                 Grid % cells_n_cells(c1)
        c2 = Grid % cells_c(i_cel, c1)
        s  = Grid % cells_f(i_cel, c1)
        fl = Flow % v_flux % n(s)

        w1 = Grid % f(s)
        if(c1.gt.c2) w1 = 1.0 - w1
        w2 = 1.0 - w1

        cfs = w1 * dens(c1) + w2 * dens(c2)
        a12 = 0.0
        a21 = 0.0

        if(c1 .lt. c2) then
          if(fl > 0.0) a21 = a21 + fl * cfs
          if(fl < 0.0) a12 = a12 - fl * cfs
          val(pos(1,s)) = val(pos(1,s)) - a12
          val(pos(2,s)) = val(pos(2,s)) - a21
        end if

        if(c1 .gt. c2) then
          if(fl > 0.0) a12 = a12 + fl * cfs
        end if

        ! Update only diaginal at c1 to avoid race conditions
        val(dia(c1)) = val(dia(c1)) + a12

      end do

    end do
    !$tf-acc loop end

  end if

  !--------------------------------------------!
  !   Viscous coefficients on the boundaries   !
  !--------------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. WALLFL  .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        a12 = visc_eff(c1) * fc(s)
        val(dia(c1)) = val(dia(c1)) + a12
      end do
      !$tf-acc loop end
     end if  ! boundary condition
  end do    ! regions

  !---------------------------------------!
  !   Upwind blending on the boundaries   !
  !---------------------------------------!
  if(Flow % u % blend_matrix) then
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c1 = Grid % faces_c(1,s)   ! inside cell
          fl = Flow % v_flux % n(s)
          val(dia(c1)) = val(dia(c1)) - min(fl, 0.0) * dens(c1)
        end do
        !$tf-acc loop end

      end if
    end do
  end if

  !-------------------------------------------------!
  !                                                 !
  !   Diagonal matrix entry for the unsteady term   !
  !                                                 !
  !-------------------------------------------------!
  if(present(dt)) then
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present, was independent
      val(dia(c)) = val(dia(c)) + dens(c) * Grid % vol(c) / dt
    end do
    !$tf-acc loop end
  end if

  !--------------------------------------------------------------!
  !   Store volume divided by central coefficient for momentum   !
  !   and refresh its buffers before discretizing the pressure   !
  !--------------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    Flow % v_m(c) = Grid % vol(c) / val(dia(c))
  end do
  !$tf-acc loop end

  ! This call is needed, the above loop goes through inside cells only
  call Grid % Exchange_Inside_Cells_Real(Flow % v_m)

  !-------------------------------------!
  !                                     !
  !   Part 1 of the under-relaxation    !
  !   (Part 2 is in Compute_Momentum)   !
  !                                     !
  !-------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$tf-acc loop end

# if T_FLOWS_DEBUG == 1
  allocate(temp(Grid % n_cells));  temp(:) = 0.0
  do c = Cells_In_Domain()  ! this is for debugging, don't do it on GPU
    ! or: temp(c) = val(dia(c))
    ! or: temp(c) = Acon % row(c+1) - Acon % row(c)
    temp(c) = flow_v_m(c)
  end do
  call Grid % Exchange_Inside_Cells_Real(temp)
  call Grid % Save_Debug_Vtu("v_m",              &
                             inside_name="v_m",  &
                             inside_cell=temp)
# endif

  call Profiler % Stop('Form_Momentum_Matrix')

  end subroutine
