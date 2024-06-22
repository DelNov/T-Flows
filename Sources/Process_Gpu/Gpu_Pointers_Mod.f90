#include "../Shared/Assert.h90"
#include "../Shared/Browse.h90"
#include "../Shared/Macros.h90"
#include "../Shared/Unused.h90"

!==============================================================================!
  module Gpu_Pointers_Mod
!------------------------------------------------------------------------------!
!   These are the links to other objects data members which are transferred    !
!   to GPUs.  Shockingly, but this is the only way to tell Nvidia Fortran      !
!   compiler that data is on the device.  This is really needed.               !
!                                                                              !
!   This module is, in essence, like an address book for Nvidia compiler.      !
!   Mind you, also a globally available address book, these data is global.    !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  implicit none
!==============================================================================!

  ! Access to Grid members
  integer, contiguous, pointer :: grid_faces_c(:,:)
  integer, contiguous, pointer :: grid_cells_n_cells(:)
  integer, contiguous, pointer :: grid_cells_c(:,:)
  integer, contiguous, pointer :: grid_cells_f(:,:)
  integer                      :: grid_n_regions,  grid_n_bnd_regions
  integer, contiguous, pointer :: grid_region_f_cell(:), grid_region_l_cell(:)
  integer, contiguous, pointer :: grid_region_f_face(:), grid_region_l_face(:)
  real,    contiguous, pointer :: grid_dx(:), grid_dy(:), grid_dz(:)
  real,    contiguous, pointer :: grid_sx(:), grid_sy(:), grid_sz(:)
  real,    contiguous, pointer :: grid_xc(:), grid_yc(:), grid_zc(:)
  real,    contiguous, pointer :: grid_d(:), grid_s(:)
  real,    contiguous, pointer :: grid_vol(:), grid_wall_dist(:)

  ! Access to Flow members
  real, contiguous, pointer :: flow_t_n(:), flow_t_o(:), flow_t_oo(:)
  real, contiguous, pointer :: flow_u_n(:), flow_u_o(:), flow_u_oo(:)
  real, contiguous, pointer :: flow_v_n(:), flow_v_o(:), flow_v_oo(:)
  real, contiguous, pointer :: flow_w_n(:), flow_w_o(:), flow_w_oo(:)
  real, contiguous, pointer :: flow_p_n(:), flow_pp_n(:), flow_v_m(:)
  real, contiguous, pointer :: flow_shear(:), flow_vort(:), flow_v_flux_n(:)
  real, contiguous, pointer :: flow_grad_c2c(:,:)
  real, contiguous, pointer :: flow_density(:),  flow_viscosity(:)
  real, contiguous, pointer :: flow_capacity(:), flow_conductivity(:)

  ! Access to Turb members
  real, contiguous, pointer :: turb_y_plus(:), turb_vis_t(:), turb_wale_v(:)

  end module
