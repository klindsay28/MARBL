module marbl_restore_mod
  !
  ! Module to generalize restoring any non-autotroph tracer
  !

  use marbl_kinds_mod, only : r8
  use marbl_kinds_mod, only : int_kind

  implicit none
  private

  !-----------------------------------------------------------------------
  !  public/private declarations
  !-----------------------------------------------------------------------

  public :: marbl_restore_compute_interior_restore
  public :: marbl_restore_compute_interior_restore_shadow

contains

!*****************************************************************************

subroutine marbl_restore_compute_interior_restore(interior_tracers, kmt,      &
                                                  interior_forcings,          &
                                                  interior_forcing_ind,       &
                                                  interior_restore)
  !
  !  restore a variable if required
  !
  use marbl_constants_mod, only : c0
  use marbl_interface_public_types, only : marbl_forcing_fields_type
  use marbl_interface_private_types, only : marbl_interior_forcing_indexing_type

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:),               intent(in) :: interior_tracers
  integer(int_kind),                           intent(in) :: kmt
  type(marbl_forcing_fields_type),             intent(in) :: interior_forcings(:)
  type(marbl_interior_forcing_indexing_type),  intent(in) :: interior_forcing_ind

  !-----------------------------------------------------------------------
  !  output variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:, :), intent(out) :: interior_restore

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------
  integer(int_kind), pointer :: restoring_inds(:)
  integer(int_kind), pointer :: inv_tau_inds(:)
  integer(int_kind) :: m, n
  !-----------------------------------------------------------------------

  interior_restore = c0
  restoring_inds => interior_forcing_ind%tracer_restore_id
  inv_tau_inds   => interior_forcing_ind%inv_tau_id

  do m=1,size(interior_forcing_ind%tracer_id)
    n = interior_forcing_ind%tracer_id(m)
    associate(restore_field => interior_forcings(restoring_inds(n))%field_1d, &
              inv_tau       =>  interior_forcings(inv_tau_inds(n))%field_1d)
      interior_restore(n,1:kmt) = (restore_field(1,1:kmt) - interior_tracers(n,1:kmt)) * inv_tau(1,1:kmt)
    end associate
  end do

end subroutine marbl_restore_compute_interior_restore

!*****************************************************************************

subroutine marbl_restore_compute_interior_restore_shadow(interior_tracers, kmt,      &
                                                         interior_forcings,          &
                                                         interior_forcing_ind,       &
                                                         marbl_tracer_indices,       &
                                                         interior_restore)
  !
  !  selectively restore shadow tracers to their non-shadow analogue
  !
  use marbl_interface_public_types, only : marbl_forcing_fields_type
  use marbl_interface_private_types, only : marbl_interior_forcing_indexing_type
  use marbl_interface_private_types, only : marbl_tracer_index_type

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:),               intent(in) :: interior_tracers
  integer(int_kind),                           intent(in) :: kmt
  type(marbl_forcing_fields_type),             intent(in) :: interior_forcings(:)
    type(marbl_interior_forcing_indexing_type),  intent(in) :: interior_forcing_ind
  type(marbl_tracer_index_type),               intent(in) :: marbl_tracer_indices

  !-----------------------------------------------------------------------
  !  input/output variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:), intent(inout) :: interior_restore

  !-----------------------------------------------------------------------

  call marbl_restore_compute_interior_restore_shadow_single( &
    interior_tracers, kmt, interior_forcings, &
    marbl_tracer_indices%po4_shadow_ind, marbl_tracer_indices%po4_ind, &
    interior_forcing_ind%normalized_POP_remin_id, &
    interior_restore)

  call marbl_restore_compute_interior_restore_shadow_single( &
    interior_tracers, kmt, interior_forcings, &
    marbl_tracer_indices%sio3_shadow_ind, marbl_tracer_indices%sio3_ind, &
    interior_forcing_ind%normalized_bSi_remin_id, &
    interior_restore)

end subroutine marbl_restore_compute_interior_restore_shadow

!*****************************************************************************

subroutine marbl_restore_compute_interior_restore_shadow_single( &
    interior_tracers, kmt, interior_forcings, &
    shadow_tracer_ind, non_shadow_tracer_ind, norm_remin_forcing_ind, &
    interior_restore)

  use marbl_interface_public_types, only : marbl_forcing_fields_type
  use marbl_settings_mod, only : parm_NK_nut_restore_invtau_peryear
  use marbl_constants_mod, only : yps

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:),               intent(in) :: interior_tracers
  integer(int_kind),                           intent(in) :: kmt
  type(marbl_forcing_fields_type),             intent(in) :: interior_forcings(:)
  integer(int_kind),                           intent(in) :: shadow_tracer_ind
  integer(int_kind),                           intent(in) :: non_shadow_tracer_ind
  integer(int_kind),                           intent(in) :: norm_remin_forcing_ind

  !-----------------------------------------------------------------------
  !  input/output variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:), intent(inout) :: interior_restore

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  integer(int_kind) :: k_surf
  real(kind=r8)     :: shadow_inv_tau
  real(kind=r8)     :: restore_tend

  !-----------------------------------------------------------------------

  k_surf = 1
  shadow_inv_tau = parm_NK_nut_restore_invtau_peryear * yps

  restore_tend = shadow_inv_tau * (interior_tracers(non_shadow_tracer_ind,k_surf) - interior_tracers(shadow_tracer_ind,k_surf))

  interior_restore(shadow_tracer_ind,k_surf) = interior_restore(shadow_tracer_ind,k_surf) + restore_tend

  interior_restore(shadow_tracer_ind,1:kmt) = interior_restore(shadow_tracer_ind,1:kmt) &
    - restore_tend * interior_forcings(norm_remin_forcing_ind)%field_1d(1,1:kmt)

end subroutine marbl_restore_compute_interior_restore_shadow_single

!*****************************************************************************

end module marbl_restore_mod
