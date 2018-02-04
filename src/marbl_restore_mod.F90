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

subroutine marbl_restore_compute_interior_restore(interior_tracers, km,       &
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
  integer(int_kind),                           intent(in) :: km
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
      interior_restore(n,:) = (restore_field(1,:) - interior_tracers(n,:)) * inv_tau(1,:)
    end associate
  end do

end subroutine marbl_restore_compute_interior_restore

!*****************************************************************************

subroutine marbl_restore_compute_interior_restore_shadow(interior_tracers, km,       &
                                                         marbl_tracer_indices,       &
                                                         interior_restore)
  !
  !  selectively restore shadow tracers to their non-shadow analogue
  !
  use marbl_constants_mod, only : yps
  use marbl_interface_private_types, only : marbl_tracer_index_type

  !-----------------------------------------------------------------------
  !  input variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:), intent(in) :: interior_tracers
  integer(int_kind),             intent(in) :: km
  type(marbl_tracer_index_type), intent(in) :: marbl_tracer_indices

  !-----------------------------------------------------------------------
  !  input/output variables
  !-----------------------------------------------------------------------

  real(kind=r8), dimension(:,:), intent(inout) :: interior_restore

  !-----------------------------------------------------------------------
  !  local variables
  !-----------------------------------------------------------------------

  integer(int_kind) :: shadow_ind, non_shadow_ind
  integer(int_kind) :: k
  real(kind=r8)     :: shadow_inv_tau

  !-----------------------------------------------------------------------
  ! SiO3, surface layer only
  !-----------------------------------------------------------------------

  shadow_ind     = marbl_tracer_indices%sio3_shadow_ind
  non_shadow_ind = marbl_tracer_indices%sio3_ind

  k = 1
  shadow_inv_tau = 10.0_r8 * yps ! 10/yr
  interior_restore(shadow_ind,k) = interior_restore(shadow_ind,k) &
    + shadow_inv_tau * (interior_tracers(non_shadow_ind,k) - interior_tracers(shadow_ind,k))

end subroutine marbl_restore_compute_interior_restore_shadow

!*****************************************************************************

end module marbl_restore_mod
