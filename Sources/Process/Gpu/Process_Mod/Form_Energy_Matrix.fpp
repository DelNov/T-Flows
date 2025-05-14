!==============================================================================!
  subroutine Form_Energy_Matrix(Process, Grid, Flow, Turb,  &
                                cond_eff, urf, dt)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)                   :: Process
  type(Grid_Type),   intent(in), target :: Grid
  type(Field_Type),              target :: Flow
  type(Turb_Type),               target :: Turb
  real                                  :: cond_eff(-Grid % n_bnd_cells &
                                                    :Grid % n_cells)
  real                                  :: urf
  real,    optional, intent(in)         :: dt       !! time step
!-----------------------------------[Locals]-----------------------------------!
  real,      contiguous, pointer :: val(:), fc(:)
  integer,   contiguous, pointer :: dia(:), pos(:,:)
  real,      contiguous, pointer :: dens_capa(:)
  integer                        :: c, s, c1, c2, i_cel, reg, nz, i
  real                           :: a12, a21, fl, cfs
# if T_FLOWS_DEBUG == 1
  real, allocatable :: temp(:)
# endif
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Form_Energy_Matrix')

  !-----------------------!
  !   Take some aliases   !
  !-----------------------!
  val => Flow % Nat % A % val
  dia => Flow % Nat % A % dia
  pos => Flow % Nat % A % pos
  fc  => Flow % Nat % A % fc
  nz  =  Flow % Nat % A % nonzeros

  call Work % Connect_Real_Cell(dens_capa)

  Assert(urf > 0.0)

  !--------------------------------------------------------------------------!
  !   Initialize density times thermal capacity and effective conductivity   !
  !--------------------------------------------------------------------------!

  !$tf-acc loop begin
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()  ! all present
    dens_capa(c) = Flow % density(c) * Flow % capacity(c)
  end do
  !$tf-acc loop end

  !-----------------------------------------------------------!
  !   Start by copying molecular viscosity to the effective   !
  !-----------------------------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()
    cond_eff(c) = Flow % conductivity(c)
  end do
  !$tf-acc loop end

  !----------------------------------------------------------------!
  !   If there is a turbulence model, add turbulent conductivity   !
  !----------------------------------------------------------------!
  if(Turb % model .ne. NO_TURBULENCE_MODEL) then
    !$tf-acc loop begin
    do c = Cells_In_Domain_And_Buffers()
      cond_eff(c) = cond_eff(c) + turb_vis_t(c) / 0.9  ! hard-coded Pr_t
    end do
    !$tf-acc loop end
  end if

  !---------------------------------------!
  !   Initialize matrix entries to zero   !
  !---------------------------------------!

  !$tf-acc loop begin
  do i = 1, nz  ! all present
    val(i) = 0.0
  end do
  !$tf-acc loop end

  !--------------------------------------------------!
  !                                                  !
  !   Compute neighbouring coefficients over cells   !
  !                                                  !
  !--------------------------------------------------!

  !------------------------------------!
  !   Coefficients inside the domain   !
  !------------------------------------!

  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)

      if(c2 .gt. 0) then

        a12 = Face_Value(s, cond_eff(c1), cond_eff(c2)) * fc(s)
        a21 = a12

        if(c1 .lt. c2) then
          val(pos(1,s)) = -a12
          val(pos(2,s)) = -a21
        end if

        ! Update only diaginal at c1 to avoid race conditions
        val(dia(c1)) = val(dia(c1)) + a12

      end if
    end do

  end do
  !$tf-acc loop end

  !---------------------------------------!
  !   Upwind blending inside the domain   !
  !---------------------------------------!
  if(Flow % t % blend_matrix) then

    !$tf-acc loop begin
    do c1 = Cells_In_Domain()  ! all present

      do i_cel = 1, Grid % cells_n_cells(c1)
        c2 = Grid % cells_c(i_cel, c1)
        s  = Grid % cells_f(i_cel, c1)
        fl = Flow % v_flux % n(s)

        if(c2 .gt. 0) then

          cfs = Face_Value(s, dens_capa(c1), dens_capa(c2))
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

        end if
      end do

    end do
    !$tf-acc loop end

  end if

  !------------------------------------!
  !   Coefficients on the boundaries   !
  !------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL    .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        a12 = cond_eff(c1) * fc(s)
        val(dia(c1)) = val(dia(c1)) + a12
      end do
      !$tf-acc loop end

    end if
  end do

  if(Flow % t % blend_matrix) then
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. INFLOW) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c1 = Grid % faces_c(1,s)   ! inside cell
          fl = Flow % v_flux % n(s)
          val(dia(c1)) = val(dia(c1)) - min(fl, 0.0) * dens_capa(c1)
        end do
        !$tf-acc loop end

      end if
    end do
  end if

  !------------------------------------!
  !                                    !
  !   Take care of the unsteady term   !
  !                                    !
  !------------------------------------!
  if(present(dt)) then
    !$tf-acc loop begin
    do c = Cells_In_Domain()  ! all present, was independent
      val(dia(c)) = val(dia(c)) + dens_capa(c) * Grid % vol(c) / dt
    end do
    !$tf-acc loop end
  end if

  !------------------------------------!
  !                                    !
  !   Part 1 of the under-relaxation   !
  !   (Part 2 is in Compute_Energy)    !
  !                                    !
  !------------------------------------!
  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present, was independent
    val(dia(c)) = val(dia(c)) / urf
  end do
  !$tf-acc loop end

  call Work % Disconnect_Real_Cell(dens_capa)

  call Profiler % Stop('Form_Energy_Matrix')

  end subroutine
