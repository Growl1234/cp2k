!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Pair-local LRI backend for k-point HFX.
!> \author SY Wang (05.2026)
! **************************************************************************************************
MODULE hfx_lri_kp
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: do_potential_short
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_lri_update_ks_kp

CONTAINS

! **************************************************************************************************
!> \brief Updates the k-point KS matrices with the pair-local LRI HFX contribution.
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_lri_update_ks_kp(qs_env, ri_data, ks_matrix, ehfx, rho_ao, &
                                   geometry_did_change, nspins, hf_fraction)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      ehfx = 0.0_dp

      IF (ri_data%hfx_pot%potential_type /= do_potential_short) THEN
         CPABORT("LRI-HFXk is currently restricted to short-range/HSE exchange")
      END IF

      ! Keep dummy arguments referenced in this scaffold.
      IF (ASSOCIATED(qs_env) .AND. ASSOCIATED(ks_matrix) .AND. ASSOCIATED(rho_ao) .AND. &
          geometry_did_change .AND. nspins > 0 .AND. hf_fraction >= 0.0_dp) THEN
         CONTINUE
      END IF

      CPABORT("LRI-HFXk backend selected, but the pair-local LRI implementation is not available yet")
   END SUBROUTINE hfx_lri_update_ks_kp

END MODULE hfx_lri_kp
