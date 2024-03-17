!==============================================================================!
  subroutine Add_Advection_Term(Proc, Flow, comp)
!------------------------------------------------------------------------------!
!   Thoroughly re-vamped for the GPU_2                                         !
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
!   Dimension of the system under consideration                                !
!     [M]{u} = {b}   [kgm/s^2]   [N]                                           !
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Proc
  type(Field_Type), target :: Flow
  integer                  :: comp
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: ui_n(:), b(:), v_flux(:), dens(:)
  real                      :: b_tmp, den_u1, den_u2, dens_f, ui_c, blend
  integer                   :: s, c1, c2, i_cel
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Proc)
!==============================================================================!

  call Profiler % Start('Add_Advection_Term')

  ! Take some aliases
  b      => Flow % Nat % b
  v_flux => Flow % v_flux
  dens   => Flow % density
  blend  =  Flow % blend

  ! Still on aliases
  if(comp .eq. 1) ui_n => Flow % u % n
  if(comp .eq. 2) ui_n => Flow % v % n
  if(comp .eq. 3) ui_n => Flow % w % n

  !-------------------------------------------!
  !   Browse through all the interior cells   !
  !      (This can be accelerted on GPU)      !
  !-------------------------------------------!

  !$acc parallel loop
  do c1 = 1, grid_n_cells - grid_n_buff_cells
    b_tmp = b(c1)
    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      if(c2 .gt. 0) then
        ui_c = 0.5 * (ui_n(c1) + ui_n(c2))  ! centered value
        ! Unit: kg / m^3 * m /s = kg / (m^2 s)
        dens_f = 0.5 * (dens(c1) + dens(c2))
        den_u1 = dens_f * ((1.0-blend) * ui_n(c1) + blend * ui_c)
        den_u2 = dens_f * ((1.0-blend) * ui_n(c2) + blend * ui_c)
        ! Unit: kg / (m^2 s) * m^3 / s = kg m / s^2 = N
        b_tmp = b_tmp - den_u1 * max(v_flux(s), 0.0) * merge(1,0, c1.lt.c2)
        b_tmp = b_tmp - den_u2 * min(v_flux(s), 0.0) * merge(1,0, c1.lt.c2)
        b_tmp = b_tmp + den_u2 * max(v_flux(s), 0.0) * merge(1,0, c1.gt.c2)
        b_tmp = b_tmp + den_u1 * min(v_flux(s), 0.0) * merge(1,0, c1.gt.c2)
      end if
    end do
    !$acc end loop
    b(c1) = b_tmp
  end do
  !$acc end parallel

  call Profiler % Stop('Add_Advection_Term')

  end subroutine
